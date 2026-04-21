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
    const bool perVariantIsoQAOnlyActive = (!forEmbeddedSim && perVariantIsoQA_ONLY);
    
    const vector<string> ueVariants = perVariantIsoQAOnlyActive
    ? vector<string>{"baseVariant"}
    : vector<string>{"noSub", "baseVariant", "variantA", "variantB"};
    const vector<string> ueLabels   = perVariantIsoQAOnlyActive
    ? vector<string>{"with UE Sub"}
    : vector<string>{"No UE Sub", "with UE Sub", "Variant A", "Variant B"};
    
    const auto& centBins = CentBins();
    if (centBins.empty())
    {
        cout << ANSI_BOLD_YEL << "[WARN] No centrality bins defined — skipping UE comparisons\n"
        << ANSI_RESET;
        return;
    }
    
    vector<CentBin> centPlotBins = centBins;
    {
        CentBin cb;
        cb.lo = 0;
        cb.hi = 20;
        cb.label = "0-20";
        cb.folder = "0_20";
        cb.suffix = "_cent_0_20";
        centPlotBins.push_back(cb);
    }
    {
        CentBin cb;
        cb.lo = 20;
        cb.hi = 60;
        cb.label = "20-60";
        cb.folder = "20_60";
        cb.suffix = "_cent_20_60";
        centPlotBins.push_back(cb);
    }
    
    auto MaterializeMergedCentralityHistogram = [&](TDirectory* dir,
                                                    const string& outName,
                                                    const vector<string>& srcNames) -> void
    {
        if (!dir) return;
        if (dir->Get(outName.c_str())) return;
        
        TH1* hMerged = nullptr;
        int nFound = 0;
        
        for (const auto& srcName : srcNames)
        {
            TH1* hSrc = dynamic_cast<TH1*>(dir->Get(srcName.c_str()));
            if (!hSrc) continue;
            ++nFound;
            
            if (!hMerged)
            {
                hMerged = CloneTH1(hSrc, outName);
                if (hMerged)
                {
                    EnsureSumw2(hMerged);
                    hMerged->Reset("ICES");
                }
            }
            if (hMerged) hMerged->Add(hSrc);
        }
        
        if (!hMerged || nFound != (int)srcNames.size())
        {
            if (hMerged) delete hMerged;
            return;
        }
        
        hMerged->SetName(outName.c_str());
        hMerged->SetTitle(outName.c_str());
        dir->Append(hMerged);
    };
    
    auto MaterializeMergedCentralityHistogramsForDir = [&](TDirectory* dir) -> void
    {
        if (!dir) return;
        
        const vector<string> histBases = {
            "h_Eiso",
            "h_Eiso_tight",
            "h_Eiso_nonTight",
            "h_EisoReco_truthSigMatched",
            "h_EisoReco_truthSigMatched_tight"
        };
        
        for (const auto& pb : PtBins())
        {
            for (const auto& base : histBases)
            {
                MaterializeMergedCentralityHistogram(
                                                     dir,
                                                     base + pb.suffix + "_cent_0_20",
                                                     {
                                                         base + pb.suffix + "_cent_0_10",
                                                         base + pb.suffix + "_cent_10_20"
                                                     }
                                                     );
                
                MaterializeMergedCentralityHistogram(
                                                     dir,
                                                     base + pb.suffix + "_cent_20_60",
                                                     {
                                                         base + pb.suffix + "_cent_20_40",
                                                         base + pb.suffix + "_cent_40_60"
                                                     }
                                                     );
            }
        }
    };
    
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
    
    // pp photonJet20 SIM file (always open, independent of sim toggles)
    TFile* fPPSim510_20 = nullptr;
    TDirectory* ppSimTop = nullptr;
    {
        const string ppSimPath = SimInputPathForSample(SimSample::kPhotonJet20);
        if (!ppSimPath.empty())
        {
            fPPSim510_20 = TFile::Open(ppSimPath.c_str(), "READ");
            if (fPPSim510_20 && !fPPSim510_20->IsZombie())
                ppSimTop = fPPSim510_20->GetDirectory(kDirSIM.c_str());
            else { if (fPPSim510_20) { fPPSim510_20->Close(); delete fPPSim510_20; fPPSim510_20 = nullptr; } }
        }
    }
    
    // pp 5+10+20 merged SIM file for standalone pythia-only iso-cut efficiency plot
    TFile* fPPSimMerged = nullptr;
    TDirectory* ppSimMergedTop = nullptr;
    {
        const string ppMergedPath = SimInputPathForSample(SimSample::kPhotonJet5And10And20Merged);
        if (!ppMergedPath.empty())
        {
            fPPSimMerged = TFile::Open(ppMergedPath.c_str(), "READ");
            if (fPPSimMerged && !fPPSimMerged->IsZombie())
                ppSimMergedTop = fPPSimMerged->GetDirectory(kDirSIM.c_str());
            else { if (fPPSimMerged) { fPPSimMerged->Close(); delete fPPSimMerged; fPPSimMerged = nullptr; } }
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
        const string layerByLayerBase = JoinPath(ueCompBase, "layerByLayerOverlays");
        const string tightNonTightBase = JoinPath(ueCompBase, "tightNonTightOverlays");
        EnsureDir(ueCompBase);
        EnsureDir(tightNonTightBase);
        if (!perVariantIsoQAOnlyActive)
        {
            EnsureDir(centralitySummaryBase);
            EnsureDir(meanIsoSummaryDir);
            EnsureDir(ptSummaryBase);
            EnsureDir(ptMeanIsoSummaryDir);
            EnsureDir(perVariantOverlayBase);
            EnsureDir(layerByLayerBase);
        }
        
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
        vector<vector<vector<double>>> gaussMean(nVarSlots, vector<vector<double>>(kNPtBins, vector<double>(centPlotBins.size(), 0.0)));
        vector<vector<vector<double>>> gaussSigma(nVarSlots, vector<vector<double>>(kNPtBins, vector<double>(centPlotBins.size(), 0.0)));
        vector<vector<vector<double>>> gaussMeanErr(nVarSlots, vector<vector<double>>(kNPtBins, vector<double>(centPlotBins.size(), 0.0)));
        vector<vector<vector<double>>> gaussSigmaErr(nVarSlots, vector<vector<double>>(kNPtBins, vector<double>(centPlotBins.size(), 0.0)));
        vector<vector<vector<bool>>>   gaussFilled(nVarSlots, vector<vector<bool>>(kNPtBins, vector<bool>(centPlotBins.size(), false)));
        
        // Restricted-range Gaussian around the mode (ALICE/ATLAS isolation approach):
        // seed from peak bin, fit asymmetrically weighting the left (signal) side,
        // truncating the right-side combinatoric tail.
        auto FitGaussianIterative = [](TH1* h, double& outMean, double& outSigma,
                                       double& outMeanErr, double& outSigmaErr) -> bool
        {
            outMean = 0.0;
            outSigma = 0.0;
            outMeanErr = 0.0;
            outSigmaErr = 0.0;
            if (!h || h->GetEntries() < 30) return false;
            
            const int maxBin = h->GetMaximumBin();
            const double peakY = h->GetBinContent(maxBin);
            if (!(peakY > 0.0)) return false;
            
            const double mode = h->GetBinCenter(maxBin);
            const double halfMax = 0.5 * peakY;
            const double hwhmToSigma = 1.0 / std::sqrt(2.0 * std::log(2.0));
            
            double hwhmL = -1.0;
            double hwhmR = -1.0;
            for (int ib = maxBin - 1; ib >= 1; --ib)
            {
                if (h->GetBinContent(ib) <= halfMax)
                {
                    hwhmL = mode - h->GetBinCenter(ib);
                    break;
                }
            }
            for (int ib = maxBin + 1; ib <= h->GetNbinsX(); ++ib)
            {
                if (h->GetBinContent(ib) <= halfMax)
                {
                    hwhmR = h->GetBinCenter(ib) - mode;
                    break;
                }
            }
            
            double sig = 0.0;
            if (hwhmL > 0.0 && hwhmR > 0.0)      sig = 0.5 * (hwhmL + hwhmR) * hwhmToSigma;
            else if (hwhmL > 0.0)                sig = hwhmL * hwhmToSigma;
            else if (hwhmR > 0.0)                sig = hwhmR * hwhmToSigma;
            else                                 sig = std::max(0.5, 0.35 * h->GetRMS());
            
            if (!(sig > 0.0) || !std::isfinite(sig)) sig = 0.5;
            
            double mu = mode;
            
            for (int pass = 0; pass < 3; ++pass)
            {
                const double nL = (pass == 0) ? 1.75 : 2.0;
                const double nR = (pass == 0) ? 1.00 : 1.25;
                const double lo = mu - nL * sig;
                const double hi = mu + nR * sig;
                
                TF1 fg(TString::Format("fg_iter_%p_%d", (void*)h, pass).Data(), "gaus", lo, hi);
                fg.SetParameters(std::max(h->GetBinContent(h->FindBin(mu)), 1.0), mu, sig);
                fg.SetParLimits(0, 0.0, 10.0 * std::max(peakY, 1.0));
                fg.SetParLimits(1, mode - 2.0 * std::max(sig, 0.5), mode + 2.0 * std::max(sig, 0.5));
                fg.SetParLimits(2, 0.10, 5.0 * std::max(sig, 0.5));
                
                const int st = h->Fit(&fg, "QNR0", "", lo, hi);
                if (st != 0 && st != 4000) continue;
                
                mu = fg.GetParameter(1);
                sig = std::abs(fg.GetParameter(2));
                if (!(sig > 0.0) || !std::isfinite(sig)) return false;
            }
            
            const double loF = mu - 2.0 * sig;
            const double hiF = mu + 1.25 * sig;
            TF1 gf(TString::Format("gf_final_%p", (void*)h).Data(), "gaus", loF, hiF);
            gf.SetParameters(std::max(h->GetBinContent(h->FindBin(mu)), 1.0), mu, sig);
            gf.SetParLimits(0, 0.0, 10.0 * std::max(peakY, 1.0));
            gf.SetParLimits(1, mode - 2.0 * std::max(sig, 0.5), mode + 2.0 * std::max(sig, 0.5));
            gf.SetParLimits(2, 0.10, 5.0 * std::max(sig, 0.5));
            
            const int st = h->Fit(&gf, "QNR0", "", loF, hiF);
            if (st != 0 && st != 4000) return false;
            
            outMean     = gf.GetParameter(1);
            outSigma    = std::abs(gf.GetParameter(2));
            outMeanErr  = gf.GetParError(1);
            outSigmaErr = gf.GetParError(2);
            
            if (!(std::isfinite(outMean) && std::isfinite(outSigma) && outSigma > 0.0)) return false;
            return true;
        };
        
        // Representative pp-only Gaussian around the mode:
        // intentionally fit a very tight core around the highest bin so the
        // plotted #mu/#sigma track the visible peak rather than the long right tail.
        auto FitGaussianRepresentativePP = [](TH1* h, double& outMean, double& outSigma,
                                              double& outMeanErr, double& outSigmaErr) -> bool
        {
            outMean = 0.0;
            outSigma = 0.0;
            outMeanErr = 0.0;
            outSigmaErr = 0.0;
            if (!h || h->GetEntries() < 10) return false;
            
            const int maxBin = h->GetMaximumBin();
            const double peakY = h->GetBinContent(maxBin);
            if (!(peakY > 0.0)) return false;
            
            const double mode = h->GetBinCenter(maxBin);
            const double binW = h->GetBinWidth(maxBin);
            
            const double leftFrac  = 0.75;
            const double rightFrac = 0.55;
            const double nLeftFrac  = std::sqrt(2.0 * std::log(1.0 / leftFrac));
            const double nRightFrac = std::sqrt(2.0 * std::log(1.0 / rightFrac));
            
            double xLeft = mode - 0.70 * binW;
            double xRight = mode + 0.80 * binW;
            bool foundLeft = false;
            bool foundRight = false;
            
            for (int ib = maxBin - 1; ib >= 1; --ib)
            {
                if (h->GetBinContent(ib) <= leftFrac * peakY)
                {
                    xLeft = h->GetBinCenter(ib);
                    foundLeft = true;
                    break;
                }
            }
            for (int ib = maxBin + 1; ib <= h->GetNbinsX(); ++ib)
            {
                if (h->GetBinContent(ib) <= rightFrac * peakY)
                {
                    xRight = h->GetBinCenter(ib);
                    foundRight = true;
                    break;
                }
            }
            
            double sigL = 0.0;
            double sigR = 0.0;
            if (foundLeft && mode > xLeft && nLeftFrac > 0.0) sigL = (mode - xLeft) / nLeftFrac;
            if (foundRight && xRight > mode && nRightFrac > 0.0) sigR = (xRight - mode) / nRightFrac;
            
            double sigSeed = 0.0;
            if (sigL > 0.0 && sigR > 0.0)      sigSeed = 0.65 * sigL + 0.35 * sigR;
            else if (sigL > 0.0)               sigSeed = sigL;
            else if (sigR > 0.0)               sigSeed = sigR;
            else                               sigSeed = std::max(0.18, 0.30 * h->GetRMS());
            
            if (!(sigSeed > 0.0) || !std::isfinite(sigSeed)) sigSeed = std::max(0.18, 0.50 * binW);
            
            double mu = mode;
            double sig = sigSeed;
            
            for (int pass = 0; pass < 3; ++pass)
            {
                const double lo = std::max(xLeft,  mu - 0.95 * sig);
                const double hi = std::min(xRight, mu + 0.85 * sig);
                if (!(hi > lo)) continue;
                
                TF1 fg(TString::Format("fg_ppasym_%p_%d", (void*)h, pass).Data(), "gaus", lo, hi);
                fg.SetParameters(std::max(h->GetBinContent(h->FindBin(mu)), 1.0), mu, sig);
                fg.SetParLimits(0, 0.0, 10.0 * std::max(peakY, 1.0));
                fg.SetParLimits(1, mode - 0.45 * std::max(sig, binW), mode + 0.40 * std::max(sig, binW));
                fg.SetParLimits(2, 0.05, 1.8 * std::max(sig, 0.12));
                
                const int st = h->Fit(&fg, "QNR0", "", lo, hi);
                if (st != 0 && st != 4000) continue;
                
                mu = fg.GetParameter(1);
                sig = std::abs(fg.GetParameter(2));
                if (!(sig > 0.0) || !std::isfinite(sig)) break;
            }
            
            const double loF = std::max(xLeft,  mu - 0.90 * sig);
            const double hiF = std::min(xRight, mu + 0.80 * sig);
            if (!(hiF > loF))
            {
                outMean     = mu;
                outSigma    = std::max(sig, 0.18);
                outMeanErr  = 0.5 * binW;
                outSigmaErr = 0.5 * binW;
                return true;
            }
            
            TF1 gf(TString::Format("gf_ppcore_%p", (void*)h).Data(), "gaus", loF, hiF);
            gf.SetParameters(std::max(h->GetBinContent(h->FindBin(mu)), 1.0), mu, sig);
            gf.SetParLimits(0, 0.0, 10.0 * std::max(peakY, 1.0));
            gf.SetParLimits(1, mode - 0.40 * std::max(sig, binW), mode + 0.35 * std::max(sig, binW));
            gf.SetParLimits(2, 0.05, 1.6 * std::max(sig, 0.12));
            
            const int st = h->Fit(&gf, "QNR0", "", loF, hiF);
            if (st == 0 || st == 4000)
            {
                outMean     = gf.GetParameter(1);
                outSigma    = std::abs(gf.GetParameter(2));
                outMeanErr  = gf.GetParError(1);
                outSigmaErr = gf.GetParError(2);
                if (std::isfinite(outMean) && std::isfinite(outSigma) && outSigma > 0.0) return true;
            }
            
            outMean     = mu;
            outSigma    = std::max(sig, 0.18);
            outMeanErr  = 0.5 * binW;
            outSigmaErr = 0.5 * binW;
            return true;
        };
        
        // Draw a Gaussian fit curve on the current pad matching a histogram's color
        auto DrawGaussFitCurve = [&](TH1* h, int color, bool usePPRepresentative = false) -> TF1*
        {
            if (!h || h->GetEntries() < (usePPRepresentative ? 10 : 30)) return nullptr;
            
            double mu, sig, muE, sigE;
            TH1* hTmp = (TH1*)h->Clone(TString::Format("hTmp_gaussFit_%p", (void*)h).Data());
            hTmp->SetDirectory(nullptr);
            
            const bool ok = usePPRepresentative
            ? FitGaussianRepresentativePP(hTmp, mu, sig, muE, sigE)
            : FitGaussianIterative(hTmp, mu, sig, muE, sigE);
            if (!ok)
            {
                delete hTmp;
                return nullptr;
            }
            
            const double fitLo = usePPRepresentative ? (mu - 1.50 * sig) : (mu - 2.0 * sig);
            const double fitHi = usePPRepresentative ? (mu + 0.75 * sig) : (mu + 1.25 * sig);
            TF1 gf(TString::Format("gf_drawSeed_%p", (void*)h).Data(), "gaus", fitLo, fitHi);
            gf.SetParameters(std::max(hTmp->GetBinContent(hTmp->FindBin(mu)), 1.0), mu, sig);
            gf.SetParLimits(0, 0.0, 10.0 * std::max(hTmp->GetBinContent(hTmp->GetMaximumBin()), 1.0));
            if (usePPRepresentative)
            {
                gf.SetParLimits(1, mu - 0.60 * std::max(sig, 0.20), mu + 0.60 * std::max(sig, 0.20));
                gf.SetParLimits(2, 0.05, 3.0 * std::max(sig, 0.20));
            }
            else
            {
                gf.SetParLimits(1, mu - 1.0 * sig, mu + 1.0 * sig);
                gf.SetParLimits(2, 0.10, 5.0 * std::max(sig, 0.5));
            }
            
            const int st = hTmp->Fit(&gf, "QNR0", "", fitLo, fitHi);
            if (st != 0 && st != 4000)
            {
                delete hTmp;
                return nullptr;
            }
            
            const double drawScale = usePPRepresentative
            ? std::max(std::abs(gf.GetParameter(2)), 0.75 * hTmp->GetBinWidth(hTmp->GetMaximumBin()))
            : std::abs(gf.GetParameter(2));
            const double drawLo = usePPRepresentative ? (gf.GetParameter(1) - 2.5 * drawScale) : fitLo;
            const double drawHi = usePPRepresentative ? (gf.GetParameter(1) + 3.5 * drawScale) : fitHi;
            
            TF1* fDraw = new TF1(TString::Format("fGauss_%p", (void*)h).Data(),
                                 "gaus", drawLo, drawHi);
            fDraw->SetParameters(gf.GetParameter(0), gf.GetParameter(1), std::abs(gf.GetParameter(2)));
            fDraw->SetLineColor(color);
            fDraw->SetLineStyle(usePPRepresentative ? 1 : 2);
            fDraw->SetLineWidth(usePPRepresentative ? 4 : 2);
            fDraw->SetNpx(usePPRepresentative ? 1200 : 400);
            fDraw->Draw("SAME");
            
            delete hTmp;
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
            if (!perVariantIsoQAOnlyActive)
            {
                EnsureDir(variantCentralitySummaryDir);
                EnsureDir(variantPtSummaryDir);
            }
            
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
            
            MaterializeMergedCentralityHistogramsForDir(aaTop);
            MaterializeMergedCentralityHistogramsForDir(incMCvarTop);
            MaterializeMergedCentralityHistogramsForDir(phoMCvarTop);
            
            // -- layer-by-layer overlays (active UE variant only, written directly under UEcomparisons/) --
            if (!perVariantIsoQAOnlyActive && !forEmbeddedSim && (H.variant == "noSub" || H.variant == "baseVariant"))
            {
                auto MergeLayerHist = [&](TDirectory* dir,
                                          const string& histBase,
                                          const vector<string>& ptSuffixes,
                                          const vector<string>& centSuffixes,
                                          const string& cloneName) -> TH1*
                {
                    if (!dir) return nullptr;
                    
                    TH1* hSum = nullptr;
                    for (const auto& ptSuf : ptSuffixes)
                    {
                        for (const auto& centSuf : centSuffixes)
                        {
                            const string hName = histBase + ptSuf + centSuf;
                            TH1* hSrc = dynamic_cast<TH1*>(dir->Get(hName.c_str()));
                            if (!hSrc) continue;
                            
                            if (!hSum)
                            {
                                hSum = CloneTH1(hSrc, cloneName.c_str());
                                if (hSum)
                                {
                                    hSum->SetDirectory(nullptr);
                                    EnsureSumw2(hSum);
                                    hSum->Reset("ICES");
                                }
                            }
                            if (hSum) hSum->Add(hSrc);
                        }
                    }
                    return hSum;
                };
                
                auto PrepareLayerHist = [&](TH1* h, bool isData, int color, int rebinFactor = 1) -> bool
                {
                    if (!h) return false;
                    
                    if (rebinFactor > 1) h->Rebin(rebinFactor);
                    EnsureSumw2(h);
                    const double integ = h->Integral(0, h->GetNbinsX() + 1);
                    if (!(integ > 0.0)) return false;
                    
                    h->Scale(1.0 / integ);
                    h->SetTitle("");
                    h->SetLineColor(color);
                    h->SetMarkerColor(color);
                    h->SetLineWidth(2);
                    h->SetFillStyle(0);
                    
                    if (isData)
                    {
                        h->SetMarkerStyle(20);
                        h->SetMarkerSize(1.0);
                    }
                    else
                    {
                        h->SetMarkerSize(0.0);
                        for (int ib = 0; ib <= h->GetNbinsX() + 1; ++ib) h->SetBinError(ib, 0.0);
                    }
                    return true;
                };
                
                const string photonLegend =
                PhotonEmbeddedShortLabel().empty()
                ? "#gamma+jet embedded"
                : TString::Format("#gamma+jet %s embedded", PhotonEmbeddedShortLabel().c_str()).Data();
                
                const string inclusiveLegend =
                InclusiveEmbeddedShortLabel().empty()
                ? "inclusive jet embedded"
                : TString::Format("inclusive jet %s embedded", InclusiveEmbeddedShortLabel().c_str()).Data();
                
                struct LayerCentCfg
                {
                    string folder;
                    vector<string> centSuffixes;
                    string label;
                };
                
                const vector<LayerCentCfg> layerCentCfgs = {
                    {"0_20",  {"_cent_0_10", "_cent_10_20"}, "0-20%"},
                    {"60_80", {"_cent_60_80"},               "60-80%"}
                };
                
                struct LayerPanelCfg
                {
                    string histBase;
                    string panelTitle;
                };
                
                const vector<LayerPanelCfg> layerPanels = {
                    {"h_Eiso_emcal",   "EMCal"},
                    {"h_Eiso_hcalin",  "IHCal"},
                    {"h_Eiso_hcalout", "OHCal"}
                };
                
                auto BuildMergedPtSuffixes = [&](int ptMin, int ptMax) -> vector<string>
                {
                    vector<string> out;
                    for (const auto& pb : PtBins())
                    {
                        if (pb.lo >= ptMin && pb.hi <= ptMax) out.push_back(pb.suffix);
                    }
                    return out;
                };
                
                struct LayerPtCfg
                {
                    string folder;
                    vector<string> ptSuffixes;
                    string label;
                };
                
                vector<LayerPtCfg> layerPtCfgs;
                layerPtCfgs.reserve(PtBins().size() + 2);
                
                for (const auto& pb : PtBins())
                {
                    LayerPtCfg cfg;
                    cfg.folder = pb.folder;
                    cfg.ptSuffixes.push_back(pb.suffix);
                    cfg.label = TString::Format("%d-%d GeV", pb.lo, pb.hi).Data();
                    layerPtCfgs.push_back(cfg);
                }
                
                {
                    LayerPtCfg cfg;
                    cfg.folder = "pT_16_35";
                    cfg.ptSuffixes = BuildMergedPtSuffixes(16, 35);
                    cfg.label = "16-35 GeV";
                    if (!cfg.ptSuffixes.empty()) layerPtCfgs.push_back(cfg);
                }
                
                {
                    LayerPtCfg cfg;
                    cfg.folder = "pT_20_35";
                    cfg.ptSuffixes = BuildMergedPtSuffixes(20, 35);
                    cfg.label = "20-35 GeV";
                    if (!cfg.ptSuffixes.empty()) layerPtCfgs.push_back(cfg);
                }
                
                for (const auto& lc : layerCentCfgs)
                {
                    const string layerVariantDir = JoinPath(layerByLayerBase, H.variant);
                    EnsureDir(layerVariantDir);
                    const string layerCentDir = JoinPath(layerVariantDir, lc.folder);
                    EnsureDir(layerCentDir);
                    
                    for (const auto& lp : layerPtCfgs)
                    {
                        TCanvas cLayer(
                                       TString::Format("c_layerByLayer_%s_%s_%s", trigAA.c_str(), lc.folder.c_str(), lp.folder.c_str()).Data(),
                                       "c_layerByLayer", 2100, 700
                                       );
                        cLayer.Divide(3, 1, 0.002, 0.002);
                        
                        vector<TH1*> keepAlive;
                        keepAlive.reserve(9);
                        vector<TObject*> keepObjs;
                        
                        for (std::size_t ipanel = 0; ipanel < layerPanels.size(); ++ipanel)
                        {
                            cLayer.cd((int)ipanel + 1);
                            gPad->SetLeftMargin(0.20);
                            gPad->SetRightMargin(0.04);
                            gPad->SetBottomMargin(0.14);
                            gPad->SetTopMargin(0.14);
                            gPad->SetTicks(1,1);
                            
                            const auto& panel = layerPanels[ipanel];
                            
                            TH1* hData = MergeLayerHist(
                                                        aaTop,
                                                        panel.histBase,
                                                        lp.ptSuffixes,
                                                        lc.centSuffixes,
                                                        TString::Format("hLayerData_%s_%s_%s_%s_%s",
                                                                        trigAA.c_str(), H.variant.c_str(),
                                                                        lc.folder.c_str(), lp.folder.c_str(),
                                                                        panel.histBase.c_str()).Data()
                                                        );
                            TH1* hPho = MergeLayerHist(
                                                       phoMCvarTop,
                                                       panel.histBase,
                                                       lp.ptSuffixes,
                                                       lc.centSuffixes,
                                                       TString::Format("hLayerPho_%s_%s_%s_%s_%s",
                                                                       trigAA.c_str(), H.variant.c_str(),
                                                                       lc.folder.c_str(), lp.folder.c_str(),
                                                                       panel.histBase.c_str()).Data()
                                                       );
                            TH1* hInc = MergeLayerHist(
                                                       incMCvarTop,
                                                       panel.histBase,
                                                       lp.ptSuffixes,
                                                       lc.centSuffixes,
                                                       TString::Format("hLayerInc_%s_%s_%s_%s_%s",
                                                                       trigAA.c_str(), H.variant.c_str(),
                                                                       lc.folder.c_str(), lp.folder.c_str(),
                                                                       panel.histBase.c_str()).Data()
                                                       );
                            
                            if (hData) keepAlive.push_back(hData);
                            if (hPho)  keepAlive.push_back(hPho);
                            if (hInc)  keepAlive.push_back(hInc);
                            
                            const int layerRebin = (ipanel == 0) ? 10 : 1;  // EMCal=10, IHCal/OHCal=1
                            const bool haveData = PrepareLayerHist(hData, true, kBlack, layerRebin);
                            const bool havePho  = PrepareLayerHist(hPho,  false, kRed + 1, layerRebin);
                            const bool haveInc  = PrepareLayerHist(hInc,  false, kBlue + 1, layerRebin);
                            
                            if (!haveData && !havePho && !haveInc)
                            {
                                TLatex tMiss;
                                tMiss.SetNDC(true);
                                tMiss.SetTextFont(42);
                                tMiss.SetTextAlign(22);
                                tMiss.SetTextSize(0.080);
                                tMiss.DrawLatex(0.50, 0.55, "MISSING");
                                
                                TLatex tPanel;
                                tPanel.SetNDC(true);
                                tPanel.SetTextFont(42);
                                tPanel.SetTextAlign(23);
                                tPanel.SetTextSize(0.060);
                                tPanel.DrawLatex(0.50, 0.96, panel.panelTitle.c_str());
                                continue;
                            }
                            
                            TH1* hFrame = havePho ? hPho : (haveInc ? hInc : hData);
                            double yMaxLayer = 0.0;
                            if (haveData) yMaxLayer = std::max(yMaxLayer, hData->GetMaximum());
                            if (havePho)  yMaxLayer = std::max(yMaxLayer, hPho->GetMaximum());
                            if (haveInc)  yMaxLayer = std::max(yMaxLayer, hInc->GetMaximum());
                            
                            hFrame->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                            hFrame->GetYaxis()->SetTitle("Normalized to Unit Area");
                            hFrame->GetXaxis()->SetTitleSize(0.055);
                            hFrame->GetYaxis()->SetTitleSize(0.055);
                            hFrame->GetXaxis()->SetLabelSize(0.045);
                            hFrame->GetYaxis()->SetLabelSize(0.045);
                            hFrame->GetYaxis()->SetTitleOffset(1.45);
                            hFrame->SetMinimum(0.0);
                            hFrame->SetMaximum((yMaxLayer > 0.0) ? (1.45 * yMaxLayer) : 1.0);
                            
                            // per-layer x-axis zoom
                            if (ipanel == 1)      hFrame->GetXaxis()->SetRangeUser(-5.0, 10.0);   // IHCal
                            else if (ipanel == 2)  hFrame->GetXaxis()->SetRangeUser(-5.0, 20.0);   // OHCal
                            
                            if (hFrame == hData) hFrame->Draw("E1");
                            else                 hFrame->Draw("hist");
                            
                            if (havePho && hPho != hFrame) hPho->Draw("hist SAME");
                            if (haveInc && hInc != hFrame) hInc->Draw("hist SAME");
                            if (haveData && hData != hFrame) hData->Draw("E1 SAME");
                            
                            // panel label (always top-right)
                            TLatex tPanel;
                            tPanel.SetNDC(true);
                            tPanel.SetTextFont(42);
                            tPanel.SetTextAlign(33);
                            tPanel.SetTextSize(0.060);
                            tPanel.DrawLatex(0.92, 0.82, panel.panelTitle.c_str());
                            
                            // sPHENIX Internal only on center pad (IHCal)
                            if (ipanel == 1)
                            {
                                TLatex tSph;
                                tSph.SetNDC(true);
                                tSph.SetTextFont(42);
                                tSph.SetTextAlign(33);
                                tSph.SetTextSize(0.042);
                                tSph.DrawLatex(0.92, 0.76, "#bf{sPHENIX} #it{Internal}");
                                tSph.SetTextSize(0.034);
                                tSph.DrawLatex(0.92, 0.71, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                            }
                            
                            // per-pad legend (top-left)
                            auto* legLayer = new TLegend(0.22, 0.72, 0.58, 0.85);
                            legLayer->SetBorderSize(0);
                            legLayer->SetFillStyle(0);
                            legLayer->SetTextFont(42);
                            legLayer->SetTextSize(0.03);
                            if (haveData) legLayer->AddEntry(hData, "Run25 Au+Au", "ep");
                            if (havePho)  legLayer->AddEntry(hPho, photonLegend.c_str(), "l");
                            if (haveInc)  legLayer->AddEntry(hInc, inclusiveLegend.c_str(), "l");
                            legLayer->Draw();
                            keepObjs.push_back(legLayer);
                            
                        }
                        
                        // -- draw canvas-wide title --
                        {
                            const string ueTag = (H.variant == "noSub") ? "no UE sub" : "with UE sub";
                            const string canvasTitle = TString::Format(
                                                                       "Run25 Au+Au vs Photon/Inclusive Jet Embedded MC, %s, p_{T}^{#gamma} = %s, %s",
                                                                       lc.label.c_str(), lp.label.c_str(), ueTag.c_str()).Data();
                            TLatex tTitle;
                            tTitle.SetNDC(true);
                            tTitle.SetTextFont(42);
                            tTitle.SetTextAlign(23);
                            tTitle.SetTextSize(0.055);
                            cLayer.cd(0);
                            tTitle.DrawLatex(0.50, 0.98, canvasTitle.c_str());
                        }
                        
                        const bool isIntegratedPt =
                        (lp.folder == "pT_16_35" || lp.folder == "pT_20_35");
                        
                        const string layerPtOutDir = isIntegratedPt
                        ? JoinPath(layerCentDir, "integratedOverPt")
                        : layerCentDir;
                        EnsureDir(layerPtOutDir);
                        
                        SaveCanvas(cLayer, JoinPath(layerPtOutDir,
                                                    TString::Format("layerByLayer_%s.png", lp.folder.c_str()).Data()));
                        
                        for (auto* hKeep : keepAlive) delete hKeep;
                        for (auto* o : keepObjs) delete o;
                    }
                }
            }
            
            // -- accumulators for <E_T^iso> vs centrality (per pT bin) --
            vector<vector<double>> vsCent_yPP(kNPtBins);
            vector<vector<double>> vsCent_eyPP(kNPtBins);
            vector<vector<double>> vsCent_yAA(kNPtBins);
            vector<vector<double>> vsCent_eyAA(kNPtBins);
            vector<vector<bool>>   vsCent_filled(kNPtBins, vector<bool>(centPlotBins.size(), false));
            vector<vector<double>> vsCent_yPPgauss(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<double>> vsCent_eyPPgauss(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<bool>>   vsCent_filledPPgauss(kNPtBins, vector<bool>(centPlotBins.size(), false));
            vector<double>         ppGaussMeanByPt(kNPtBins, 0.0);
            vector<double>         ppGaussMeanErrByPt(kNPtBins, 0.0);
            vector<double>         ppGaussSigmaByPt(kNPtBins, 0.0);
            vector<double>         ppGaussSigmaErrByPt(kNPtBins, 0.0);
            vector<bool>           ppGaussFilledByPt(kNPtBins, false);
            for (int ip = 0; ip < kNPtBins; ++ip)
            {
                vsCent_yPP[ip].resize(centPlotBins.size(), 0.0);
                vsCent_eyPP[ip].resize(centPlotBins.size(), 0.0);
                vsCent_yAA[ip].resize(centPlotBins.size(), 0.0);
                vsCent_eyAA[ip].resize(centPlotBins.size(), 0.0);
            }
            
            if (ppTop)
            {
                for (int ipt = 0; ipt < kNPtBins; ++ipt)
                {
                    const PtBin& b = PtBins()[ipt];
                    const string hPPName = "h_Eiso" + b.suffix;
                    TH1* hPPsrc = dynamic_cast<TH1*>(ppTop->Get(hPPName.c_str()));
                    if (!hPPsrc) continue;
                    
                    double ppMu = 0.0, ppMuE = 0.0, ppSig = 0.0, ppSigE = 0.0;
                    TH1* hPPFit = CloneTH1(hPPsrc,
                                           TString::Format("hPP_meanVsCent_fit_%s_%s_%zu",
                                                           trigAA.c_str(), H.variant.c_str(), (std::size_t)ipt).Data());
                    bool ppGaussOk = false;
                    if (hPPFit)
                    {
                        EnsureSumw2(hPPFit);
                        ppGaussOk = FitGaussianIterative(hPPFit, ppMu, ppSig, ppMuE, ppSigE);
                        delete hPPFit;
                    }
                    
                    if (!ppGaussOk) continue;
                    
                    ppGaussMeanByPt[ipt]     = ppMu;
                    ppGaussMeanErrByPt[ipt]  = ppMuE;
                    ppGaussSigmaByPt[ipt]    = ppSig;
                    ppGaussSigmaErrByPt[ipt] = ppSigE;
                    ppGaussFilledByPt[ipt]   = true;
                }
            }
            
            // MC mean iso accumulators [ipt][ic] for vs-centrality and vs-pT summary plots
            vector<vector<double>> vsCent_yIncMC(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<double>> vsCent_eyIncMC(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<double>> vsCent_yPhoMC(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<double>> vsCent_eyPhoMC(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<bool>>   vsCent_filledIncMC(kNPtBins, vector<bool>(centPlotBins.size(), false));
            vector<vector<bool>>   vsCent_filledPhoMC(kNPtBins, vector<bool>(centPlotBins.size(), false));
            
            // Signal (tight data + truth-sig-matched MC) mean iso accumulators
            vector<vector<double>> vsCent_yTightData(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<double>> vsCent_eyTightData(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<bool>>   vsCent_filledTightData(kNPtBins, vector<bool>(centPlotBins.size(), false));
            vector<vector<double>> vsCent_ySigIncMC(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<double>> vsCent_eySigIncMC(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<bool>>   vsCent_filledSigIncMC(kNPtBins, vector<bool>(centPlotBins.size(), false));
            vector<vector<double>> vsCent_ySigPhoMC(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<double>> vsCent_eySigPhoMC(kNPtBins, vector<double>(centPlotBins.size(), 0.0));
            vector<vector<bool>>   vsCent_filledSigPhoMC(kNPtBins, vector<bool>(centPlotBins.size(), false));
            
            for (std::size_t ic = 0; ic < centPlotBins.size(); ++ic)
            {
                const auto& cb = centPlotBins[ic];
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
                    
                    // Gaussian fit on raw pp histogram:
                    // use the exact same pp Gaussian already used for
                    // centralitySummaryPerPt/baseVariant/ppIso_counts_*
                    if (ppGaussFilledByPt[ipt])
                    {
                        vsCent_yPPgauss[ipt][ic]      = ppGaussMeanByPt[ipt];
                        vsCent_eyPPgauss[ipt][ic]     = ppGaussMeanErrByPt[ipt];
                        vsCent_filledPPgauss[ipt][ic] = true;
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
                                    hTd->Scale(1.0, "width");
                                    hNTd->Scale(1.0, "width");
                                    hSig->Scale(1.0, "width");
                                    hTd->Rebin(10);
                                    hNTd->Rebin(10);
                                    hSig->Rebin(10);
                                    {
                                        const double meanIsoET = hTd->GetMean();
                                        const double rmsIsoET  = hTd->GetRMS();
                                        const double normCut   = meanIsoET + 1.75 * rmsIsoET;
                                        
                                        const int binNormLo = hTd->FindBin(normCut);
                                        const int binNormHi = hTd->GetNbinsX();
                                        
                                        const double normTight    = hTd->Integral(binNormLo, binNormHi);
                                        const double normNonTight = hNTd->Integral(binNormLo, binNormHi);
                                        if (normNonTight > 0.0) hNTd->Scale(normTight / normNonTight);
                                        
                                        const double dataSignalYield = hTd->Integral() - hNTd->Integral();
                                        const double mcSignalYield   = hSig->Integral();
                                        if (mcSignalYield > 0.0) hSig->Scale(dataSignalYield / mcSignalYield);
                                        
                                        hSig->SetLineColor(kBlue-9); hSig->SetFillColorAlpha(kBlue-9,0.5); hSig->SetFillStyle(1001); hSig->SetLineWidth(1); hSig->SetMarkerSize(0.0);
                                        hNTd->SetLineColor(kPink-4); hNTd->SetFillColorAlpha(kPink-4,0.5); hNTd->SetFillStyle(1001); hNTd->SetLineWidth(1); hNTd->SetMarkerSize(0.0);
                                        for (int ib=0; ib<=hSig->GetNbinsX()+1; ++ib) hSig->SetBinError(ib,0.0);
                                        for (int ib=0; ib<=hNTd->GetNbinsX()+1; ++ib) hNTd->SetBinError(ib,0.0);
                                        hTd->SetLineColor(kBlack); hTd->SetMarkerColor(kBlack); hTd->SetMarkerStyle(20); hTd->SetMarkerSize(1.0); hTd->SetLineWidth(2); hTd->SetFillStyle(0);
                                        TH1* hSigTopPP = (TH1*)hSig->Clone(TString::Format("hSigTop_ppSig_%s_%s", cb.folder.c_str(), b.folder.c_str()).Data());
                                        hSigTopPP->SetDirectory(nullptr); hSigTopPP->Add(hNTd);
                                        const double ym = std::max(hTd->GetMaximum(), hSigTopPP->GetMaximum());
                                        TCanvas cSB(TString::Format("c_ppSigBkg_%s_%s", cb.folder.c_str(), b.folder.c_str()).Data(), "c_ppSigBkg", 900, 700);
                                        ApplyCanvasMargins1D(cSB); cSB.cd();
                                        hSigTopPP->SetTitle(""); hSigTopPP->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                                        hSigTopPP->GetYaxis()->SetTitle("Counts / Bin Width");
                                        hTd->GetXaxis()->SetTitleSize(0.055); hTd->GetYaxis()->SetTitleSize(0.055);
                                        hTd->GetXaxis()->SetLabelSize(0.045); hTd->GetYaxis()->SetLabelSize(0.045);
                                        hTd->GetYaxis()->SetTitleOffset(1.30);
                                        hSigTopPP->GetXaxis()->SetTitleSize(0.055); hSigTopPP->GetYaxis()->SetTitleSize(0.055);
                                        hSigTopPP->GetXaxis()->SetLabelSize(0.045); hSigTopPP->GetYaxis()->SetLabelSize(0.045);
                                        hSigTopPP->GetYaxis()->SetTitleOffset(1.30);
                                        hSigTopPP->SetMinimum(0.0); hSigTopPP->SetMaximum((ym>0)?(1.3*ym):1.0);
                                        hSigTopPP->SetLineColor(kBlue-9); hSigTopPP->SetFillColorAlpha(kBlue-9,0.5); hSigTopPP->SetFillStyle(1001);
                                        hSigTopPP->Draw("hist");
                                        hNTd->Draw("hist SAME");
                                        hTd->Draw("E1 SAME");
                                        TLatex tSBtitle; tSBtitle.SetNDC(true); tSBtitle.SetTextFont(42); tSBtitle.SetTextAlign(23); tSBtitle.SetTextSize(0.040);
                                        tSBtitle.DrawLatex(0.50, 0.97,
                                                           "E_{T}^{iso}, p+p data versus Photon+Jet (5+10+20) Pythia");
                                        TLegend lg(0.56,0.65,0.92,0.88); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextFont(42); lg.SetTextSize(0.032);
                                        lg.AddEntry(hTd, "Data (Tight)", "ep"); lg.AddEntry(hNTd, "Data (Non-tight)", "f"); lg.AddEntry(hSig, "Signal Embedded MC", "f"); lg.Draw();
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
                                        delete hSigTopPP;
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
                                const string hSigTightName = "h_EisoReco_truthSigMatched_tight" + b.suffix + cb.suffix;
                                TH1* hTdSrc  = dynamic_cast<TH1*>(aaTop->Get(hTName.c_str()));
                                TH1* hNTdSrc = dynamic_cast<TH1*>(aaTop->Get(hNTName.c_str()));
                                TH1* hSigSrc = dynamic_cast<TH1*>(mcCfg.mcTop->Get(hSigName.c_str()));
                                TH1* hSigTightSrc = dynamic_cast<TH1*>(mcCfg.mcTop->Get(hSigTightName.c_str()));
                                TH1* hRecoTightSrc = dynamic_cast<TH1*>(mcCfg.mcTop->Get(hTName.c_str()));
                                if (hTdSrc && hNTdSrc && (hSigSrc || hSigTightSrc || hRecoTightSrc))
                                {
                                    auto DrawSignalBkgOverlay = [&](TH1* hSigTemplateSrc,
                                                                    const string& cloneTag,
                                                                    const string& outName,
                                                                    const string& sigLegend) -> void
                                    {
                                        if (!hSigTemplateSrc) return;
                                        
                                        TH1* hTd  = CloneTH1(hTdSrc,  TString::Format("hTd_aaSig_%s_%s_%s_%s_%s", cloneTag.c_str(), mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                                        TH1* hNTd = CloneTH1(hNTdSrc, TString::Format("hNTd_aaSig_%s_%s_%s_%s_%s", cloneTag.c_str(), mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                                        TH1* hSig = CloneTH1(hSigTemplateSrc, TString::Format("hSig_aaSig_%s_%s_%s_%s_%s", cloneTag.c_str(), mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                                        if (!(hTd && hNTd && hSig))
                                        {
                                            if (hTd) delete hTd;
                                            if (hNTd) delete hNTd;
                                            if (hSig) delete hSig;
                                            return;
                                        }
                                        
                                        hTd->Scale(1.0, "width");
                                        hNTd->Scale(1.0, "width");
                                        hSig->Scale(1.0, "width");
                                        hTd->Rebin(10);
                                        hNTd->Rebin(10);
                                        hSig->Rebin(10);
                                        
                                        const double meanIsoET = hTd->GetMean();
                                        const double rmsIsoET  = hTd->GetRMS();
                                        const double normCut   = meanIsoET + 1.75 * rmsIsoET;
                                        
                                        const int binNormLo = hTd->FindBin(normCut);
                                        const int binNormHi = hTd->GetNbinsX();
                                        
                                        const double normTight    = hTd->Integral(binNormLo, binNormHi);
                                        const double normNonTight = hNTd->Integral(binNormLo, binNormHi);
                                        if (normNonTight > 0.0) hNTd->Scale(normTight / normNonTight);
                                        
                                        const double dataSignalYield = hTd->Integral() - hNTd->Integral();
                                        const double mcSignalYield   = hSig->Integral();
                                        if (mcSignalYield > 0.0) hSig->Scale(dataSignalYield / mcSignalYield);
                                        
                                        hSig->SetLineColor(kBlue-9); hSig->SetFillColorAlpha(kBlue-9,0.5); hSig->SetFillStyle(1001); hSig->SetLineWidth(1); hSig->SetMarkerSize(0.0);
                                        hNTd->SetLineColor(kPink-4); hNTd->SetFillColorAlpha(kPink-4,0.5); hNTd->SetFillStyle(1001); hNTd->SetLineWidth(1); hNTd->SetMarkerSize(0.0);
                                        for (int ib=0; ib<=hSig->GetNbinsX()+1; ++ib) hSig->SetBinError(ib,0.0);
                                        for (int ib=0; ib<=hNTd->GetNbinsX()+1; ++ib) hNTd->SetBinError(ib,0.0);
                                        hTd->SetLineColor(kBlack); hTd->SetMarkerColor(kBlack); hTd->SetMarkerStyle(20); hTd->SetMarkerSize(1.0); hTd->SetLineWidth(2); hTd->SetFillStyle(0);
                                        TH1* hSigTop = (TH1*)hSig->Clone(TString::Format("hSigTop_aaSig_%s_%s_%s_%s_%s", cloneTag.c_str(), mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                                        hSigTop->SetDirectory(nullptr);
                                        hSigTop->Add(hNTd);
                                        const double ym = std::max(hTd->GetMaximum(), hSigTop->GetMaximum());
                                        TCanvas cSB(TString::Format("c_aaSigBkg_%s_%s_%s_%s_%s", cloneTag.c_str(), mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(), "c_aaSigBkg", 900, 700);
                                        ApplyCanvasMargins1D(cSB); cSB.cd();
                                        hSigTop->SetTitle(""); hSigTop->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                                        hSigTop->GetYaxis()->SetTitle("Counts / Bin Width");
                                        hSigTop->GetXaxis()->SetTitleSize(0.055); hSigTop->GetYaxis()->SetTitleSize(0.055);
                                        hSigTop->GetXaxis()->SetLabelSize(0.045); hSigTop->GetYaxis()->SetLabelSize(0.045);
                                        hSigTop->GetYaxis()->SetTitleOffset(1.30);
                                        hSigTop->SetMinimum(0.0); hSigTop->SetMaximum((ym>0)?(1.3*ym):1.0);
                                        hSigTop->SetLineColor(kBlue-9); hSigTop->SetFillColorAlpha(kBlue-9,0.5); hSigTop->SetFillStyle(1001);
                                        hSigTop->Draw("hist");
                                        hNTd->Draw("hist SAME");
                                        hTd->Draw("E1 SAME");
                                        TLatex tSBtitle; tSBtitle.SetNDC(true); tSBtitle.SetTextFont(42); tSBtitle.SetTextAlign(23); tSBtitle.SetTextSize(0.040);
                                        tSBtitle.DrawLatex(0.50, 0.97,
                                                           TString::Format("E_{T}^{iso} overlay: AuAu data vs %s Embedded MC", mcCfg.titleTag.c_str()).Data());
                                        TLegend lg(0.56,0.65,0.92,0.88); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextFont(42); lg.SetTextSize(0.032);
                                        lg.AddEntry(hTd, "Data (Tight)", "ep"); lg.AddEntry(hNTd, "Data (Non-tight)", "f"); lg.AddEntry(hSig, sigLegend.c_str(), "f"); lg.Draw();
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
                                        SaveCanvas(cSB, JoinPath(sigOvDir, outName));
                                        
                                        delete hSigTop;
                                        delete hTd;
                                        delete hNTd;
                                        delete hSig;
                                    };
                                    
                                    DrawSignalBkgOverlay(hSigSrc,
                                                         "inclusiveTruthMatched",
                                                         "Eiso_sigBkg_overlay.png",
                                                         "Signal Embedded MC");
                                    DrawSignalBkgOverlay(hSigTightSrc,
                                                         "tightTaggedTruthMatched",
                                                         "Eiso_sigBkg_overlay_TIGHT_TAGGED.png",
                                                         "Signal Embedded MC (Tight)");
                                    DrawSignalBkgOverlay(hRecoTightSrc,
                                                         "recoTightMC",
                                                         "Eiso_sigBkg_overlay_RECO_TIGHT_MC.png",
                                                         "Reco tight MC");
                                }
                            }
                        }
                    }
                    // ── Merged pT 20–35 bin: combine last N bins where lo >= 20 ──                {
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
                                    TH1* hTdMSrc  = MergeHists(aaTop,      "h_Eiso_tight",    cb.suffix, "hTd_merged_"+mcCfg.folder+"_"+cb.folder);
                                    TH1* hNTdMSrc = MergeHists(aaTop,      "h_Eiso_nonTight", cb.suffix, "hNTd_merged_"+mcCfg.folder+"_"+cb.folder);
                                    TH1* hSigMSrc = MergeHists(mcCfg.mcTop,"h_EisoReco_truthSigMatched", cb.suffix, "hSig_merged_"+mcCfg.folder+"_"+cb.folder);
                                    TH1* hSigTightMSrc = MergeHists(mcCfg.mcTop,"h_EisoReco_truthSigMatched_tight", cb.suffix, "hSigTight_merged_"+mcCfg.folder+"_"+cb.folder);
                                    TH1* hRecoTightMSrc = MergeHists(mcCfg.mcTop,"h_Eiso_tight", cb.suffix, "hRecoTight_merged_"+mcCfg.folder+"_"+cb.folder);
                                    if (hTdMSrc && hNTdMSrc && (hSigMSrc || hSigTightMSrc || hRecoTightMSrc))
                                    {
                                        auto DrawMergedSignalBkgOverlay = [&](TH1* hSigTemplateSrc,
                                                                              const string& cloneTag,
                                                                              const string& outName,
                                                                              const string& sigLegend) -> void
                                        {
                                            if (!hSigTemplateSrc) return;
                                            
                                            TH1* hTdM  = CloneTH1(hTdMSrc,  TString::Format("hTdM_aaSig_%s_%s_%s", cloneTag.c_str(), mcCfg.folder.c_str(), cb.folder.c_str()).Data());
                                            TH1* hNTdM = CloneTH1(hNTdMSrc, TString::Format("hNTdM_aaSig_%s_%s_%s", cloneTag.c_str(), mcCfg.folder.c_str(), cb.folder.c_str()).Data());
                                            TH1* hSigM = CloneTH1(hSigTemplateSrc, TString::Format("hSigM_aaSig_%s_%s_%s", cloneTag.c_str(), mcCfg.folder.c_str(), cb.folder.c_str()).Data());
                                            if (!(hTdM && hNTdM && hSigM))
                                            {
                                                if (hTdM) delete hTdM;
                                                if (hNTdM) delete hNTdM;
                                                if (hSigM) delete hSigM;
                                                return;
                                            }
                                            
                                            hTdM->Scale(1.0, "width");
                                            hNTdM->Scale(1.0, "width");
                                            hSigM->Scale(1.0, "width");
                                            hTdM->Rebin(10);
                                            hNTdM->Rebin(10);
                                            hSigM->Rebin(10);
                                            
                                            const double meanIsoETM = hTdM->GetMean();
                                            const double rmsIsoETM  = hTdM->GetRMS();
                                            const double normCutM   = meanIsoETM + 1.75 * rmsIsoETM;
                                            
                                            const int binNormLoM = hTdM->FindBin(normCutM);
                                            const int binNormHiM = hTdM->GetNbinsX();
                                            
                                            const double normTightM    = hTdM->Integral(binNormLoM, binNormHiM);
                                            const double normNonTightM = hNTdM->Integral(binNormLoM, binNormHiM);
                                            if (normNonTightM > 0.0) hNTdM->Scale(normTightM / normNonTightM);
                                            
                                            const double dataSignalYieldM = hTdM->Integral() - hNTdM->Integral();
                                            const double mcSignalYieldM   = hSigM->Integral();
                                            if (mcSignalYieldM > 0.0) hSigM->Scale(dataSignalYieldM / mcSignalYieldM);
                                            
                                            hSigM->SetLineColor(kBlue-9); hSigM->SetFillColorAlpha(kBlue-9,0.5); hSigM->SetFillStyle(1001); hSigM->SetLineWidth(1); hSigM->SetMarkerSize(0.0);
                                            hNTdM->SetLineColor(kPink-4); hNTdM->SetFillColorAlpha(kPink-4,0.5); hNTdM->SetFillStyle(1001); hNTdM->SetLineWidth(1); hNTdM->SetMarkerSize(0.0);
                                            for (int ib=0;ib<=hSigM->GetNbinsX()+1;++ib) hSigM->SetBinError(ib,0.0);
                                            for (int ib=0;ib<=hNTdM->GetNbinsX()+1;++ib) hNTdM->SetBinError(ib,0.0);
                                            hTdM->SetLineColor(kBlack); hTdM->SetMarkerColor(kBlack); hTdM->SetMarkerStyle(20); hTdM->SetMarkerSize(1.0); hTdM->SetLineWidth(2); hTdM->SetFillStyle(0);
                                            TH1* hSigTopM = (TH1*)hSigM->Clone(TString::Format("hSigTopM_aaSig_%s_%s_%s", cloneTag.c_str(), mcCfg.folder.c_str(), cb.folder.c_str()).Data());
                                            hSigTopM->SetDirectory(nullptr); hSigTopM->Add(hNTdM);
                                            const double ymM = std::max(hTdM->GetMaximum(),hSigTopM->GetMaximum());
                                            TCanvas cSBm(TString::Format("c_merged_sigBkg_%s_%s_%s", cloneTag.c_str(), mcCfg.folder.c_str(), cb.folder.c_str()).Data(), "c_merged_sigBkg",900,700); ApplyCanvasMargins1D(cSBm); cSBm.cd();
                                            hSigTopM->SetTitle(""); hSigTopM->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]"); hSigTopM->GetYaxis()->SetTitle("Counts / Bin Width");
                                            hSigTopM->GetXaxis()->SetTitleSize(0.055); hSigTopM->GetYaxis()->SetTitleSize(0.055); hSigTopM->GetXaxis()->SetLabelSize(0.045); hSigTopM->GetYaxis()->SetLabelSize(0.045); hSigTopM->GetYaxis()->SetTitleOffset(1.30);
                                            hSigTopM->SetMinimum(0.0); hSigTopM->SetMaximum((ymM>0)?(1.3*ymM):1.0);
                                            hSigTopM->SetLineColor(kBlue-9); hSigTopM->SetFillColorAlpha(kBlue-9,0.5); hSigTopM->SetFillStyle(1001);
                                            hSigTopM->Draw("hist"); hNTdM->Draw("hist SAME"); hTdM->Draw("E1 SAME");
                                            TLatex tSBtM; tSBtM.SetNDC(true); tSBtM.SetTextFont(42); tSBtM.SetTextAlign(23); tSBtM.SetTextSize(0.040);
                                            tSBtM.DrawLatex(0.50,0.97,TString::Format("E_{T}^{iso} overlay: AuAu data vs %s Embedded MC",mcCfg.titleTag.c_str()).Data());
                                            TLegend lgSBm(0.56,0.65,0.92,0.88); lgSBm.SetBorderSize(0); lgSBm.SetFillStyle(0); lgSBm.SetTextFont(42); lgSBm.SetTextSize(0.032);
                                            lgSBm.AddEntry(hTdM,"Data (Tight)","ep"); lgSBm.AddEntry(hNTdM,"Data (Non-tight)","f"); lgSBm.AddEntry(hSigM,sigLegend.c_str(),"f"); lgSBm.Draw();
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
                                            SaveCanvas(cSBm,JoinPath(sigOvDirM,outName));
                                            delete hSigTopM;
                                            delete hTdM;
                                            delete hNTdM;
                                            delete hSigM;
                                        };
                                        
                                        DrawMergedSignalBkgOverlay(hSigMSrc,
                                                                   "inclusiveTruthMatched",
                                                                   "Eiso_sigBkg_overlay.png",
                                                                   "Signal Embedded MC");
                                        DrawMergedSignalBkgOverlay(hSigTightMSrc,
                                                                   "tightTaggedTruthMatched",
                                                                   "Eiso_sigBkg_overlay_TIGHT_TAGGED.png",
                                                                   "Signal Embedded MC (Tight)");
                                        DrawMergedSignalBkgOverlay(hRecoTightMSrc,
                                                                   "recoTightMC",
                                                                   "Eiso_sigBkg_overlay_RECO_TIGHT_MC.png",
                                                                   "Reco tight MC");
                                    }
                                    if (hTdMSrc) delete hTdMSrc;
                                    if (hNTdMSrc) delete hNTdMSrc;
                                    if (hSigMSrc) delete hSigMSrc;
                                    if (hSigTightMSrc) delete hSigTightMSrc;
                                    if (hRecoTightMSrc) delete hRecoTightMSrc;
                                }
                            }
                        }
                    }
                }
                {
                    vector<double> xPtGauss;
                    vector<double> exPtGauss;
                    vector<double> yPPGauss;
                    vector<double> eyPPGauss;
                    vector<double> yAAGauss;
                    vector<double> eyAAGauss;
                    
                    for (int ipt = 0; ipt < kNPtBins; ++ipt)
                    {
                        if (!ppGaussFilledByPt[ipt] || !gaussFilled[ivH][ipt][ic]) continue;
                        
                        xPtGauss.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                        exPtGauss.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                        yPPGauss.push_back(ppGaussMeanByPt[ipt]);
                        eyPPGauss.push_back(ppGaussMeanErrByPt[ipt]);
                        yAAGauss.push_back(gaussMean[ivH][ipt][ic]);
                        eyAAGauss.push_back(gaussMeanErr[ivH][ipt][ic]);
                    }
                    
                    if (!xPtGauss.empty())
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
                        for (std::size_t i = 0; i < xPtGauss.size(); ++i)
                        {
                            yMin = std::min(yMin, std::min(yPPGauss[i] - eyPPGauss[i], yAAGauss[i] - eyAAGauss[i]));
                            yMax = std::max(yMax, std::max(yPPGauss[i] + eyPPGauss[i], yAAGauss[i] + eyAAGauss[i]));
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
                        hFrame.GetYaxis()->SetTitle("#mu^{Gauss} [GeV]");
                        hFrame.GetXaxis()->SetTitleSize(0.055);
                        hFrame.GetYaxis()->SetTitleSize(0.055);
                        hFrame.GetXaxis()->SetLabelSize(0.045);
                        hFrame.GetYaxis()->SetLabelSize(0.045);
                        hFrame.GetYaxis()->SetTitleOffset(1.15);
                        hFrame.Draw();
                        
                        TGraphErrors gPP((int)xPtGauss.size(), &xPtGauss[0], &yPPGauss[0], &exPtGauss[0], &eyPPGauss[0]);
                        gPP.SetLineWidth(2);
                        gPP.SetLineColor(kRed + 1);
                        gPP.SetMarkerColor(kRed + 1);
                        gPP.SetMarkerStyle(24);
                        gPP.SetMarkerSize(1.2);
                        gPP.Draw("PE1 SAME");
                        
                        TGraphErrors gAA((int)xPtGauss.size(), &xPtGauss[0], &yAAGauss[0], &exPtGauss[0], &eyAAGauss[0]);
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
                                         TString::Format("#mu^{Gauss} overlay pp and %d-%d%% Cent AuAu", cb.lo, cb.hi).Data());
                        
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
                                lgS.AddEntry(&gTDs, "Data (Tight)", "ep");
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
                        if (!vsCent_filledPPgauss[ipt][ic] || !gaussFilled[ivH][ipt][ic]) continue;
                        const auto& cb = centBins[ic];
                        xCent.push_back(0.5 * (cb.lo + cb.hi));
                        exCent.push_back(0.5 * (cb.hi - cb.lo));
                        yCentPP.push_back(vsCent_yPPgauss[ipt][ic]);
                        eyCentPP.push_back(vsCent_eyPPgauss[ipt][ic]);
                        yCentAA.push_back(gaussMean[ivH][ipt][ic]);
                        eyCentAA.push_back(gaussMeanErr[ivH][ipt][ic]);
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
                    hFrame.GetYaxis()->SetTitle("#mu^{Gauss} [GeV]");
                    hFrame.GetXaxis()->SetTitleSize(0.055);
                    hFrame.GetYaxis()->SetTitleSize(0.055);
                    hFrame.GetXaxis()->SetLabelSize(0.045);
                    hFrame.GetYaxis()->SetLabelSize(0.045);
                    hFrame.GetYaxis()->SetTitleOffset(1.15);
                    hFrame.Draw();
                    
                    double meanPP = 0.0;
                    for (double v : yCentPP) meanPP += v;
                    if (!yCentPP.empty()) meanPP /= (double)yCentPP.size();
                    TLine lPP(centLo, meanPP, centHi, meanPP);
                    lPP.SetLineColor(kBlack);
                    lPP.SetLineStyle(2);
                    lPP.SetLineWidth(2);
                    lPP.Draw("SAME");
                    
                    TGraphErrors gAA((int)xCent.size(), &xCent[0], &yCentAA[0], &exCent[0], &eyCentAA[0]);
                    gAA.SetLineWidth(2);
                    gAA.SetLineColor(kBlack);
                    gAA.SetMarkerColor(kBlack);
                    gAA.SetMarkerStyle(20);
                    gAA.SetMarkerSize(1.2);
                    gAA.Draw("PE1 SAME");
                    
                    TLegend leg(0.58, 0.17, 0.92, 0.32);
                    leg.SetBorderSize(0);
                    leg.SetFillStyle(0);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.032);
                    leg.AddEntry(&lPP, "pp data", "l");
                    leg.AddEntry(&gAA, TString::Format("AuAu data (%s)", H.label.c_str()).Data(), "ep");
                    leg.Draw();
                    
                    TLatex tTitle;
                    tTitle.SetNDC(true);
                    tTitle.SetTextFont(42);
                    tTitle.SetTextAlign(23);
                    tTitle.SetTextSize(0.045);
                    tTitle.DrawLatex(0.50, 0.965,
                                     TString::Format("#mu^{Gauss} vs Centrality  p_{T}^{#gamma} %d-%d GeV", b.lo, b.hi).Data());
                    
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextAlign(13);
                    t.SetTextSize(0.034);
                    t.DrawLatex(0.19, 0.88, trigDisplayLabel.c_str());
                    t.DrawLatex(0.19, 0.84, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                    t.DrawLatex(0.19, 0.80, H.label.c_str());
                    
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
                                    lgS.AddEntry(&gTDs, "Data (Tight)", "ep");
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
                        { double integ = hAA->Integral(); if (integ > 0.0) hAA->Scale(1.0 / integ); }
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
                    
                    bool havePeripheralGauss = false;
                    double peripheralMu = 0.0, peripheralMuErr = 0.0;
                    double peripheralSigma = 0.0, peripheralSigmaErr = 0.0;
                    int peripheralColor = kOrange + 1;
                    for (std::size_t icSub = 0; icSub < centBins.size(); ++icSub)
                    {
                        const auto& cbSub = centBins[icSub];
                        if (cbSub.lo == 60 && cbSub.hi == 80 && gaussFilled[ivH][ipt][icSub])
                        {
                            havePeripheralGauss = true;
                            peripheralMu       = gaussMean[ivH][ipt][icSub];
                            peripheralMuErr    = gaussMeanErr[ivH][ipt][icSub];
                            peripheralSigma    = gaussSigma[ivH][ipt][icSub];
                            peripheralSigmaErr = gaussSigmaErr[ivH][ipt][icSub];
                            peripheralColor    = (icSub < 8) ? centOvColors[icSub] : (kOrange + 1);
                            break;
                        }
                    }
                    
                    bool havePPGauss = false;
                    double ppMuOverlay = 0.0, ppMuErrOverlay = 0.0;
                    double ppSigmaOverlay = 0.0, ppSigmaErrOverlay = 0.0;
                    TH1* hPPOverlayFit = nullptr;
                    if (!forEmbeddedSim && ppTop)
                    {
                        const string hPPName = "h_Eiso" + b.suffix;
                        TH1* hPPSrc = dynamic_cast<TH1*>(ppTop->Get(hPPName.c_str()));
                        if (hPPSrc)
                        {
                            hPPOverlayFit = CloneTH1(
                                                     hPPSrc,
                                                     TString::Format("hPP_overlayFit_%s_%s_%s",
                                                                     trigAA.c_str(), H.variant.c_str(), b.folder.c_str()).Data()
                                                     );
                            if (hPPOverlayFit)
                            {
                                hPPOverlayFit->Rebin(10);
                                EnsureSumw2(hPPOverlayFit);
                                havePPGauss = FitGaussianRepresentativePP(
                                                                          hPPOverlayFit,
                                                                          ppMuOverlay, ppSigmaOverlay,
                                                                          ppMuErrOverlay, ppSigmaErrOverlay
                                                                          );
                            }
                        }
                    }
                    
                    // --- Collect Gaussian means & sigmas per centrality for subpanels ---
                    vector<double> subX, subEX, subY, subEY;
                    vector<double> subSigY, subSigEY;
                    for (std::size_t icSub = 0; icSub < centBins.size(); ++icSub)
                    {
                        if (!gaussFilled[ivH][ipt][icSub]) continue;
                        const auto& cbSub = centBins[icSub];
                        subX.push_back(0.5 * (cbSub.lo + cbSub.hi));
                        subEX.push_back(0.5 * (cbSub.hi - cbSub.lo));
                        subY.push_back(gaussMean[ivH][ipt][icSub]);
                        subEY.push_back(gaussMeanErr[ivH][ipt][icSub]);
                        subSigY.push_back(gaussSigma[ivH][ipt][icSub]);
                        subSigEY.push_back(gaussSigmaErr[ivH][ipt][icSub]);
                    }
                    
                    auto DrawCentralityOverlayTop = [&](vector<TF1*>& outCurves)
                    {
                        hCents[0]->SetTitle("");
                        hCents[0]->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                        hCents[0]->GetXaxis()->SetTitleSize(0.048);
                        hCents[0]->GetXaxis()->SetLabelSize(0.04);
                        hCents[0]->GetXaxis()->SetTitleOffset(0.75);
                        hCents[0]->GetYaxis()->SetTitle("Normalized to Unit Area");
                        
                        hCents[0]->GetYaxis()->SetTitleSize(0.050);
                        hCents[0]->GetYaxis()->SetLabelSize(0.050);
                        hCents[0]->GetYaxis()->SetTitleOffset(1.15);
                        hCents[0]->GetXaxis()->SetRangeUser(-10.0, 50.0);
                        hCents[0]->SetMinimum(0.0);
                        {
                            const double yScale = (H.variant == "noSub") ? 1.20
                            : (H.variant == "baseVariant") ? 1.15
                            : 1.25;
                            hCents[0]->SetMaximum((yMaxCO > 0.0) ? (yScale * yMaxCO) : 1.0);
                        }
                        hCents[0]->Draw("E1");
                        for (std::size_t ih = 1; ih < hCents.size(); ++ih) hCents[ih]->Draw("E1 SAME");
                        
                        outCurves.clear();
                        for (std::size_t ih = 0; ih < hCents.size(); ++ih)
                        {
                            TF1* fc = DrawGaussFitCurve(hCents[ih], hCents[ih]->GetLineColor());
                            if (fc) outCurves.push_back(fc);
                        }
                        
                        const double yLineMax = hCents[0]->GetMaximum();
                        if (havePeripheralGauss && H.variant == "baseVariant")
                        {
                            TLine lPer(peripheralMu, 0.0, peripheralMu, yLineMax);
                            lPer.SetLineColor(peripheralColor);
                            lPer.SetLineStyle(2);
                            lPer.SetLineWidth(2);
                            lPer.DrawClone();
                        }
                        
                        const bool legMidRight = (H.variant == "noSub");
                        TLegend legCO(legMidRight ? 0.55 : 0.60,
                                      legMidRight ? 0.57 : 0.40,
                                      legMidRight ? 0.86 : 0.90,
                                      legMidRight ? 0.75 : 0.62);
                        legCO.SetColumnSeparation(-0.05);
                        legCO.SetBorderSize(0);
                        legCO.SetFillStyle(0);
                        legCO.SetTextFont(42);
                        legCO.SetTextSize(legMidRight ? 0.034 : 0.042);
                        legCO.SetNColumns(2);
                        for (std::size_t ih = 0; ih < hCents.size(); ++ih)
                            legCO.AddEntry(hCents[ih], cLabels[ih].c_str(), "ep");
                        legCO.DrawClone();
                        
                        TLatex tTitleCO;
                        tTitleCO.SetNDC(true);
                        tTitleCO.SetTextFont(42);
                        tTitleCO.SetTextAlign(23);
                        tTitleCO.SetTextSize(0.052);
                        tTitleCO.DrawLatex(0.50, 0.99,
                                           TString::Format("Au+Au centrality overlays, p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                        
                        TLatex tInfoCO;
                        tInfoCO.SetNDC(true);
                        tInfoCO.SetTextFont(42);
                        tInfoCO.SetTextAlign(33);
                        tInfoCO.SetTextSize((H.variant == "noSub") ? 0.036 : 0.040);
                        tInfoCO.DrawLatex(0.90, 0.88,
                                          forEmbeddedSim
                                          ? (embeddedMode == 2 ? activeEmbeddedSimFolder.substr(0, activeEmbeddedSimFolder.find("_SIM")).c_str()
                                             : SimSampleLabel(CurrentSimSample()).c_str())
                                          : trigDisplayLabel.c_str());
                        tInfoCO.DrawLatex(0.90, 0.84, H.label.c_str());
                        {
                            const double coneRValCO = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                            tInfoCO.DrawLatex(0.90, 0.80, TString::Format("#DeltaR_{cone} < %.1f", coneRValCO).Data());
                        }
                        {
                            TLatex tSph;
                            tSph.SetNDC(true);
                            tSph.SetTextFont(42);
                            tSph.SetTextAlign(33);
                            tSph.SetTextSize(0.052);
                            const double sphY1  = (H.variant == "noSub") ? 0.56 : 0.39;
                            const double sphX1  = (H.variant == "noSub") ? 0.93 : 0.90;
                            tSph.DrawLatex(sphX1, sphY1,        "#bf{sPHENIX} #it{Internal}");
                            tSph.SetTextSize(0.042);
                            tSph.DrawLatex(sphX1, sphY1 - 0.06, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                        }
                        
                        if (havePeripheralGauss && H.variant == "baseVariant")
                        {
                            TLatex tPer;
                            tPer.SetNDC(true);
                            tPer.SetTextFont(42);
                            tPer.SetTextAlign(33);
                            tPer.SetTextSize(0.032);
                            tPer.SetTextColor(peripheralColor);
                            tPer.DrawLatex(0.90, 0.74,
                                           TString::Format("#mu^{Gauss} (60-80%%) = %.3f #pm %.4f GeV", peripheralMu, peripheralMuErr).Data());
                            tPer.DrawLatex(0.90, 0.7,
                                           TString::Format("#sigma^{Gauss} (60-80%%) = %.3f #pm %.4f GeV", peripheralSigma, peripheralSigmaErr).Data());
                        }
                    };
                    
                    TCanvas cCO(
                                TString::Format("c_isoCentOverlay_%s_%s_%s",
                                                trigAA.c_str(), H.variant.c_str(), b.folder.c_str()).Data(),
                                "c_isoCentOverlay", 900, 1000
                                );
                    cCO.cd();
                    
                    // Upper pad: histogram overlay
                    TPad* padUp = new TPad(TString::Format("padUp_%s_%s", H.variant.c_str(), b.folder.c_str()).Data(),
                                           "padUp", 0.0, 0.36, 1.0, 1.0);
                    padUp->SetBottomMargin(0.10);
                    padUp->SetLeftMargin(0.14);
                    padUp->SetRightMargin(0.04);
                    padUp->SetTopMargin(0.08);
                    padUp->Draw();
                    padUp->cd();
                    
                    vector<TF1*> gaussCurvesC;
                    DrawCentralityOverlayTop(gaussCurvesC);
                    
                    // Middle pad: Gaussian mean vs centrality
                    cCO.cd();
                    TPad* padLo = new TPad(TString::Format("padLo_%s_%s", H.variant.c_str(), b.folder.c_str()).Data(),
                                           "padLo", 0.0, 0.21, 1.0, 0.36);
                    padLo->SetTopMargin(0.02);
                    padLo->SetBottomMargin(0.00);
                    padLo->SetLeftMargin(0.14);
                    padLo->SetRightMargin(0.04);
                    padLo->Draw();
                    padLo->cd();
                    
                    TH1F* hFrSub = nullptr;
                    TGraphErrors* gSub = nullptr;
                    if (!subX.empty())
                    {
                        double yLoSub = 1e30, yHiSub = -1e30;
                        for (std::size_t is = 0; is < subY.size(); ++is)
                        {
                            yLoSub = std::min(yLoSub, subY[is] - subEY[is]);
                            yHiSub = std::max(yHiSub, subY[is] + subEY[is]);
                        }
                        const double subPadFrac = (H.variant == "noSub") ? 0.15 : 0.25;
                        const double subPad = (yHiSub > yLoSub) ? subPadFrac * (yHiSub - yLoSub) : 1.0;
                        const double subMin = (H.variant == "noSub") ? 0.0 : (yLoSub - subPad);
                        const double subMax = (H.variant == "noSub") ? (yHiSub + 0.10 * (yHiSub - yLoSub)) : (yHiSub + subPad);
                        
                        hFrSub = new TH1F(TString::Format("hFrSub_centOv_%s_%s", H.variant.c_str(), b.folder.c_str()).Data(),
                                          "", 100, centBins.front().lo, centBins.back().hi);
                        hFrSub->SetDirectory(nullptr);
                        hFrSub->SetStats(0);
                        hFrSub->SetMinimum(subMin);
                        hFrSub->SetMaximum(subMax);
                        hFrSub->GetYaxis()->SetTitle("#mu^{Gauss}[GeV]");
                        hFrSub->GetYaxis()->SetTitleSize(0.19);
                        hFrSub->GetYaxis()->SetTitleOffset(0.28);
                        hFrSub->GetYaxis()->SetLabelSize(0.14);
                        hFrSub->GetYaxis()->SetNdivisions(505);
                        hFrSub->GetXaxis()->SetTitle("");
                        hFrSub->GetXaxis()->SetTitleSize(0.0);
                        hFrSub->GetXaxis()->SetLabelSize(0.0);
                        hFrSub->GetXaxis()->SetTickLength(0.0);
                        hFrSub->Draw();
                        
                        vector<double> exZeroSub(subX.size(), 0.0);
                        if (H.variant != "noSub")
                        {
                            TLine lineSub(centBins.front().lo, 1.0, centBins.back().hi, 1.0);
                            lineSub.SetLineColor(kBlack);
                            lineSub.SetLineStyle(2);
                            lineSub.SetLineWidth(1);
                            lineSub.DrawClone();
                        }
                        
                        gSub = new TGraphErrors((int)subX.size(), &subX[0], &subY[0], &exZeroSub[0], &subEY[0]);
                        gSub->SetMarkerStyle(20);
                        gSub->SetMarkerSize(1.0);
                        gSub->SetMarkerColor(kBlack);
                        gSub->SetLineColor(kBlack);
                        gSub->SetLineWidth(2);
                        gSub->Draw("PE1 SAME");
                    }
                    
                    // Bottom pad: Gaussian sigma vs centrality
                    cCO.cd();
                    TPad* padBot = new TPad(TString::Format("padBot_%s_%s", H.variant.c_str(), b.folder.c_str()).Data(),
                                            "padBot", 0.0, 0.0, 1.0, 0.21);
                    padBot->SetTopMargin(0.00);
                    padBot->SetBottomMargin(0.30);
                    padBot->SetLeftMargin(0.14);
                    padBot->SetRightMargin(0.04);
                    padBot->Draw();
                    padBot->cd();
                    
                    TH1F* hFrSig = nullptr;
                    TGraphErrors* gSig = nullptr;
                    if (!subX.empty())
                    {
                        double yLoSig = 1e30, yHiSig = -1e30;
                        for (std::size_t is = 0; is < subSigY.size(); ++is)
                        {
                            yLoSig = std::min(yLoSig, subSigY[is] - subSigEY[is]);
                            yHiSig = std::max(yHiSig, subSigY[is] + subSigEY[is]);
                        }
                        const double sigPadFrac = (H.variant == "noSub") ? 0.15 : 0.25;
                        const double sigPad = (yHiSig > yLoSig) ? sigPadFrac * (yHiSig - yLoSig) : 1.0;
                        const double sigMin = (H.variant == "noSub") ? 0.0 : (yLoSig - sigPad);
                        const double sigMax = (H.variant == "noSub") ? (yHiSig + 0.10 * (yHiSig - yLoSig)) : (yHiSig + sigPad);
                        
                        hFrSig = new TH1F(TString::Format("hFrSig_centOv_%s_%s", H.variant.c_str(), b.folder.c_str()).Data(),
                                          "", 100, centBins.front().lo, centBins.back().hi);
                        hFrSig->SetDirectory(nullptr);
                        hFrSig->SetStats(0);
                        hFrSig->SetMinimum(sigMin);
                        hFrSig->SetMaximum(sigMax);
                        hFrSig->GetYaxis()->SetTitle("#sigma^{Gauss}[GeV]");
                        hFrSig->GetYaxis()->SetTitleSize(0.14);
                        hFrSig->GetYaxis()->SetTitleOffset(0.32);
                        hFrSig->GetYaxis()->SetLabelSize(0.10);
                        hFrSig->GetYaxis()->SetNdivisions(505);
                        hFrSig->GetXaxis()->SetTitle("Centrality [%]");
                        hFrSig->GetXaxis()->SetTitleSize(0.12);
                        hFrSig->GetXaxis()->SetTitleOffset(0.95);
                        hFrSig->GetXaxis()->SetLabelSize(0.10);
                        hFrSig->Draw();
                        
                        vector<double> exZeroSig(subX.size(), 0.0);
                        gSig = new TGraphErrors((int)subX.size(), &subX[0], &subSigY[0], &exZeroSig[0], &subSigEY[0]);
                        gSig->SetMarkerStyle(20);
                        gSig->SetMarkerSize(1.0);
                        gSig->SetMarkerColor(kBlack);
                        gSig->SetLineColor(kBlack);
                        gSig->SetLineWidth(2);
                        gSig->Draw("PE1 SAME");
                    }
                    
                    cCO.Modified();
                    cCO.Update();
                    if (!perVariantIsoQAOnlyActive)
                        SaveCanvas(cCO, JoinPath(variantCentralitySummaryDir,
                                                 TString::Format("isoCentOverlay_counts_%s.png", b.folder.c_str()).Data()));
                    
                    TCanvas cCOTopOnly(
                                       TString::Format("c_isoCentOverlayTopOnly_%s_%s_%s",
                                                       trigAA.c_str(), H.variant.c_str(), b.folder.c_str()).Data(),
                                       "c_isoCentOverlayTopOnly", 900, 700
                                       );
                    ApplyCanvasMargins1D(cCOTopOnly);
                    cCOTopOnly.cd();
                    
                    vector<TF1*> gaussCurvesTopOnly;
                    DrawCentralityOverlayTop(gaussCurvesTopOnly);
                    if (H.variant == "baseVariant" && b.lo == 12 && b.hi == 14)
                        hCents[0]->GetXaxis()->SetRangeUser(-10.0, 20.0);
                    cCOTopOnly.Modified();
                    cCOTopOnly.Update();
                    if (!perVariantIsoQAOnlyActive)
                        SaveCanvas(cCOTopOnly, JoinPath(variantCentralitySummaryDir,
                                                        TString::Format("isoCentOverlay_counts_topOnly_%s.png", b.folder.c_str()).Data()));
                    for (auto* f : gaussCurvesTopOnly) delete f;
                    
                    if (gSig) delete gSig;
                    if (hFrSig) delete hFrSig;
                    if (gSub) delete gSub;
                    if (hFrSub) delete hFrSub;
                    
                    for (auto* f : gaussCurvesC) delete f;
                    for (auto* h : hCents) delete h;
                }
            }
            
            // ── Standalone pp E_T^iso per pT bin (baseVariant only) ──
            if (!forEmbeddedSim && H.variant == "baseVariant" && ppTop)
            {
                for (int ipt = 0; ipt < kNPtBins; ++ipt)
                {
                    const PtBin& b = PtBins()[ipt];
                    const string hPPName = "h_Eiso" + b.suffix;
                    TH1* hPPsrc = dynamic_cast<TH1*>(ppTop->Get(hPPName.c_str()));
                    if (!hPPsrc) continue;
                    
                    TH1* hPP = CloneTH1(hPPsrc,
                                        TString::Format("hPP_standalone_%s_%s",
                                                        trigAA.c_str(), b.folder.c_str()).Data());
                    if (!hPP) continue;
                    
                    EnsureSumw2(hPP);
                    hPP->SetLineColor(kBlack);
                    hPP->SetMarkerColor(kBlack);
                    hPP->SetMarkerStyle(20);
                    hPP->SetMarkerSize(0.9);
                    hPP->SetLineWidth(2);
                    hPP->SetFillStyle(0);
                    
                    // Gaussian fit: use the exact same fitter as Au+Au, but do NOT rebin for pp
                    double ppMu = 0.0, ppMuE = 0.0, ppSig = 0.0, ppSigE = 0.0;
                    TH1* hPPFit = CloneTH1(hPPsrc,
                                           TString::Format("hPP_standalone_fit_%s_%s",
                                                           trigAA.c_str(), b.folder.c_str()).Data());
                    bool ppGaussOk = false;
                    if (hPPFit)
                    {
                        EnsureSumw2(hPPFit);
                        ppGaussOk = FitGaussianIterative(hPPFit, ppMu, ppSig, ppMuE, ppSigE);
                    }
                    
                    TCanvas cPP(
                                TString::Format("c_ppIso_%s_%s", trigAA.c_str(), b.folder.c_str()).Data(),
                                "c_ppIso", 900, 700);
                    cPP.cd();
                    gPad->SetLeftMargin(0.14);
                    gPad->SetRightMargin(0.04);
                    gPad->SetBottomMargin(0.14);
                    gPad->SetTopMargin(0.08);
                    gPad->SetTicks(1,1);
                    
                    hPP->SetTitle("");
                    hPP->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                    hPP->GetYaxis()->SetTitle("Counts");
                    hPP->GetXaxis()->SetTitleSize(0.048);
                    hPP->GetYaxis()->SetTitleSize(0.050);
                    hPP->GetXaxis()->SetLabelSize(0.040);
                    hPP->GetYaxis()->SetLabelSize(0.040);
                    hPP->GetYaxis()->SetTitleOffset(1.15);
                    hPP->GetXaxis()->SetRangeUser(-10.0, 20.0);
                    hPP->SetMinimum(0.0);
                    hPP->SetMaximum(1.15 * hPP->GetMaximum());
                    hPP->Draw("E1");
                    
                    TF1* fDraw = nullptr;
                    if (ppGaussOk && hPPFit)
                    {
                        fDraw = DrawGaussFitCurve(hPPFit, kBlack, false);
                        if (fDraw)
                        {
                            const double fitStart = fDraw->GetXmin();
                            const double fitEnd   = fDraw->GetXmax();
                            const double newFitEnd = ppMu + 0.85 * (fitEnd - ppMu);
                            fDraw->SetRange(fitStart, newFitEnd);
                        }
                        hPP->Draw("E1 SAME");
                    }
                    
                    TLatex tTitle;
                    tTitle.SetNDC(true);
                    tTitle.SetTextFont(42);
                    tTitle.SetTextAlign(23);
                    tTitle.SetTextSize(0.048);
                    tTitle.DrawLatex(0.50, 0.99,
                                     TString::Format("p+p isolation, p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                    
                    TLatex tSph;
                    tSph.SetNDC(true);
                    tSph.SetTextFont(42);
                    tSph.SetTextAlign(33);
                    tSph.SetTextSize(0.042);
                    tSph.DrawLatex(0.88, 0.56, "#bf{sPHENIX} #it{Internal}");
                    tSph.SetTextSize(0.036);
                    tSph.DrawLatex(0.88, 0.51, "p+p  #sqrt{s} = 200 GeV");
                    
                    if (ppGaussOk)
                    {
                        TLine lPP(ppMu, 0.0, ppMu, hPP->GetMaximum());
                        lPP.SetLineColor(kBlack);
                        lPP.SetLineStyle(2);
                        lPP.SetLineWidth(2);
                        lPP.DrawClone();
                        
                        TLatex tGauss;
                        tGauss.SetNDC(true);
                        tGauss.SetTextFont(42);
                        tGauss.SetTextAlign(33);
                        tGauss.SetTextSize(0.034);
                        tGauss.DrawLatex(0.88, 0.44,
                                         TString::Format("#mu^{Gauss} = %.3f #pm %.4f GeV", ppMu, ppMuE).Data());
                        tGauss.DrawLatex(0.88, 0.39,
                                         TString::Format("#sigma^{Gauss} = %.3f #pm %.4f GeV", ppSig, ppSigE).Data());
                    }
                    
                    {
                        const double coneRVal = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                        TLatex tInfo;
                        tInfo.SetNDC(true);
                        tInfo.SetTextFont(42);
                        tInfo.SetTextAlign(33);
                        tInfo.SetTextSize(0.036);
                        tInfo.DrawLatex(0.88, 0.88, "Photon 4 GeV + MBD NS #geq 1");
                        tInfo.DrawLatex(0.88, 0.83, TString::Format("#DeltaR_{cone} < %.1f", coneRVal).Data());
                    }
                    
                    if (!perVariantIsoQAOnlyActive)
                        SaveCanvas(cPP, JoinPath(variantCentralitySummaryDir,
                                                 TString::Format("ppIso_counts_%s.png", b.folder.c_str()).Data()));
                    
                    if (fDraw) delete fDraw;
                    if (hPPFit) delete hPPFit;
                    delete hPP;
                }
            }
            
            // ── pT overlay of E_T^iso (unnormalized counts) per centrality bin, this variant ──
            
            {
                const int ptOvColors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1,
                    kViolet+1, kYellow+2, kSpring+2, kAzure+1,
                    kTeal+1, kPink+1};
                for (std::size_t ic = 0; ic < centPlotBins.size(); ++ic)
                {
                    const auto& cb = centPlotBins[ic];
                    const string centSubDir = JoinPath(variantDir, cb.folder);
                    EnsureDir(centSubDir);
                    
                    vector<TH1*> hPts;
                    vector<string> ptLabels;
                    double yMaxPO = 0.0;
                    
                    for (int ipt = 0; ipt < kNPtBins; ++ipt)
                    {
                        const PtBin& b = PtBins()[ipt];
                        if (b.lo < 10) continue;   // skip 8-10 GeV bin
                        const string hAAName = "h_Eiso" + b.suffix + cb.suffix;
                        TH1* hAAsrc = dynamic_cast<TH1*>(aaTop->Get(hAAName.c_str()));
                        if (!hAAsrc) continue;
                        
                        TH1* hAA = CloneTH1(hAAsrc,
                                            TString::Format("hAA_ptOvlay_%s_%s_%s_%s",
                                                            trigAA.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                        if (!hAA) continue;
                        
                        hAA->Rebin(10);
                        EnsureSumw2(hAA);
                        { double integ = hAA->Integral(); if (integ > 0.0) hAA->Scale(1.0 / integ); }
                        const int ci = ((int)hPts.size() < 11) ? ptOvColors[hPts.size()] : kBlack;
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
                    
                    // --- Collect Gaussian means & sigmas per p_{T} bin for subpanels ---
                    vector<double> subPtX, subPtEX, subPtY, subPtEY;
                    vector<double> subPtSigY, subPtSigEY;
                    for (int ipt2 = 0; ipt2 < kNPtBins; ++ipt2)
                    {
                        if (kPtEdges[(std::size_t)ipt2] < 10) continue;   // skip 8-10 GeV bin
                        if (!gaussFilled[ivH][ipt2][ic]) continue;
                        subPtX.push_back(0.5 * (kPtEdges[(std::size_t)ipt2] + kPtEdges[(std::size_t)ipt2 + 1]));
                        subPtEX.push_back(0.5 * (kPtEdges[(std::size_t)ipt2 + 1] - kPtEdges[(std::size_t)ipt2]));
                        subPtY.push_back(gaussMean[ivH][ipt2][ic]);
                        subPtEY.push_back(gaussMeanErr[ivH][ipt2][ic]);
                        subPtSigY.push_back(gaussSigma[ivH][ipt2][ic]);
                        subPtSigEY.push_back(gaussSigmaErr[ivH][ipt2][ic]);
                    }
                    
                    TCanvas cPO(
                                TString::Format("c_isoPtOverlay_%s_%s_%s",
                                                trigAA.c_str(), H.variant.c_str(), cb.folder.c_str()).Data(),
                                "c_isoPtOverlay", 900, 1000
                                );
                    cPO.cd();
                    
                    // Upper pad: histogram overlay
                    TPad* padUpPO = new TPad(TString::Format("padUpPO_%s_%s", H.variant.c_str(), cb.folder.c_str()).Data(),
                                             "padUpPO", 0.0, 0.36, 1.0, 1.0);
                    padUpPO->SetBottomMargin(0.10);
                    padUpPO->SetLeftMargin(0.14);
                    padUpPO->SetRightMargin(0.04);
                    padUpPO->SetTopMargin(0.08);
                    padUpPO->Draw();
                    padUpPO->cd();
                    
                    hPts[0]->SetTitle("");
                    hPts[0]->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                    hPts[0]->GetXaxis()->SetTitleSize(0.048);
                    hPts[0]->GetXaxis()->SetLabelSize(0.04);
                    hPts[0]->GetXaxis()->SetTitleOffset(0.75);
                    hPts[0]->GetYaxis()->SetTitle("Normalized to Unit Area");
                    hPts[0]->GetYaxis()->SetTitleSize(0.060);
                    hPts[0]->GetYaxis()->SetLabelSize(0.050);
                    hPts[0]->GetYaxis()->SetTitleOffset(1.15);
                    hPts[0]->GetXaxis()->SetRangeUser(-10.0, 50.0);
                    hPts[0]->SetMinimum(0.0);
                    {
                        const double yScale = (H.variant == "noSub") ? 1.35
                        : (H.variant == "baseVariant") ? 1.15
                        : 1.25;
                        hPts[0]->SetMaximum((yMaxPO > 0.0) ? (yScale * yMaxPO) : 1.0);
                    }
                    hPts[0]->Draw("E1");
                    for (std::size_t ih = 1; ih < hPts.size(); ++ih) hPts[ih]->Draw("E1 SAME");
                    
                    // Gaussian fit curves per pT bin
                    vector<TF1*> gaussCurvesP;
                    for (std::size_t ih = 0; ih < hPts.size(); ++ih)
                    {
                        TF1* fp = DrawGaussFitCurve(hPts[ih], hPts[ih]->GetLineColor());
                        if (fp) gaussCurvesP.push_back(fp);
                    }
                    
                    TLegend legPO(0.48, 0.45, 0.92, 0.73);
                    legPO.SetColumnSeparation(-0.05);
                    legPO.SetBorderSize(0);
                    legPO.SetFillStyle(0);
                    legPO.SetTextFont(42);
                    legPO.SetTextSize(0.036);
                    legPO.SetNColumns(3);
                    for (std::size_t ih = 0; ih < hPts.size(); ++ih)
                        legPO.AddEntry(hPts[ih], ptLabels[ih].c_str(), "ep");
                    legPO.Draw();
                    
                    TLatex tTitlePO;
                    tTitlePO.SetNDC(true);
                    tTitlePO.SetTextFont(42);
                    tTitlePO.SetTextAlign(23);
                    tTitlePO.SetTextSize(0.052);
                    tTitlePO.DrawLatex(0.50, 0.99,
                                       TString::Format("Au+Au (counts), p_{T}^{#gamma} overlays, cent = %d-%d%%", cb.lo, cb.hi).Data());
                    
                    TLatex tInfoPO;
                    tInfoPO.SetNDC(true);
                    tInfoPO.SetTextFont(42);
                    tInfoPO.SetTextAlign(33);
                    tInfoPO.SetTextSize(0.040);
                    tInfoPO.DrawLatex(0.90, 0.88,
                                      forEmbeddedSim
                                      ? (embeddedMode == 2 ? activeEmbeddedSimFolder.substr(0, activeEmbeddedSimFolder.find("_SIM")).c_str()
                                         : SimSampleLabel(CurrentSimSample()).c_str())
                                      : trigDisplayLabel.c_str());
                    tInfoPO.DrawLatex(0.90, 0.84, H.label.c_str());
                    {
                        const double coneRValPO = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                        tInfoPO.DrawLatex(0.90, 0.80, TString::Format("#DeltaR_{cone} < %.1f", coneRValPO).Data());
                    }
                    {
                        TLatex tSphPO;
                        tSphPO.SetNDC(true);
                        tSphPO.SetTextFont(42);
                        tSphPO.SetTextAlign(33);
                        tSphPO.SetTextSize(0.052);
                        const double sphY1PO = (H.variant == "noSub") ? 0.58 : 0.34;
                        tSphPO.DrawLatex(0.90, sphY1PO,        "#bf{sPHENIX} #it{Internal}");
                        tSphPO.SetTextSize(0.042);
                        tSphPO.DrawLatex(0.90, sphY1PO - 0.06, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                    }
                    
                    const double subPtLoEdge = kPtEdges[1];
                    const double subPtHiEdge = kPtEdges[(std::size_t)kNPtBins];
                    
                    // Middle pad: Gaussian mean vs pT
                    cPO.cd();
                    TPad* padMidPO = new TPad(TString::Format("padMidPO_%s_%s", H.variant.c_str(), cb.folder.c_str()).Data(),
                                              "padMidPO", 0.0, 0.21, 1.0, 0.36);
                    padMidPO->SetTopMargin(0.02);
                    padMidPO->SetBottomMargin(0.00);
                    padMidPO->SetLeftMargin(0.14);
                    padMidPO->SetRightMargin(0.04);
                    padMidPO->Draw();
                    padMidPO->cd();
                    
                    TH1F* hFrSubPtMean = nullptr;
                    TGraphErrors* gSubPtMean = nullptr;
                    if (!subPtX.empty())
                    {
                        double yLoPt = 1e30, yHiPt = -1e30;
                        for (std::size_t is = 0; is < subPtY.size(); ++is)
                        {
                            yLoPt = std::min(yLoPt, subPtY[is] - subPtEY[is]);
                            yHiPt = std::max(yHiPt, subPtY[is] + subPtEY[is]);
                        }
                        const double subPadPt = (yHiPt > yLoPt) ? 0.25 * (yHiPt - yLoPt) : 1.0;
                        
                        hFrSubPtMean = new TH1F(TString::Format("hFrSubMean_ptOv_%s_%s", H.variant.c_str(), cb.folder.c_str()).Data(),
                                                "", 100, subPtLoEdge, subPtHiEdge);
                        hFrSubPtMean->SetDirectory(nullptr);
                        hFrSubPtMean->SetStats(0);
                        hFrSubPtMean->SetMinimum(yLoPt - subPadPt);
                        hFrSubPtMean->SetMaximum(yHiPt + subPadPt);
                        hFrSubPtMean->GetYaxis()->SetTitle("#mu^{Gauss}[GeV]");
                        hFrSubPtMean->GetYaxis()->SetTitleSize(0.19);
                        hFrSubPtMean->GetYaxis()->SetTitleOffset(0.28);
                        hFrSubPtMean->GetYaxis()->SetLabelSize(0.12);
                        hFrSubPtMean->GetYaxis()->SetLabelOffset(0.008);
                        hFrSubPtMean->GetYaxis()->SetNdivisions(505);
                        hFrSubPtMean->GetXaxis()->SetNdivisions(505);
                        hFrSubPtMean->GetXaxis()->SetTitle("");
                        hFrSubPtMean->GetXaxis()->SetTitleSize(0.0);
                        hFrSubPtMean->GetXaxis()->SetLabelSize(0.0);
                        hFrSubPtMean->GetXaxis()->SetTickLength(0.0);
                        hFrSubPtMean->Draw();
                        
                        TLine lineSubPt(subPtLoEdge, 1.0, subPtHiEdge, 1.0);
                        lineSubPt.SetLineColor(kBlack);
                        lineSubPt.SetLineStyle(2);
                        lineSubPt.SetLineWidth(1);
                        lineSubPt.DrawClone();
                        
                        vector<double> exZeroSubPt(subPtX.size(), 0.0);
                        gSubPtMean = new TGraphErrors((int)subPtX.size(), &subPtX[0], &subPtY[0], &exZeroSubPt[0], &subPtEY[0]);
                        gSubPtMean->SetMarkerStyle(20);
                        gSubPtMean->SetMarkerSize(1.0);
                        gSubPtMean->SetMarkerColor(kBlack);
                        gSubPtMean->SetLineColor(kBlack);
                        gSubPtMean->SetLineWidth(2);
                        gSubPtMean->Draw("PE1 SAME");
                    }
                    
                    // Bottom pad: Gaussian sigma vs pT
                    cPO.cd();
                    TPad* padBotPO = new TPad(TString::Format("padBotPO_%s_%s", H.variant.c_str(), cb.folder.c_str()).Data(),
                                              "padBotPO", 0.0, 0.0, 1.0, 0.21);
                    padBotPO->SetTopMargin(0.00);
                    padBotPO->SetBottomMargin(0.30);
                    padBotPO->SetLeftMargin(0.14);
                    padBotPO->SetRightMargin(0.04);
                    padBotPO->Draw();
                    padBotPO->cd();
                    
                    TH1F* hFrSubPtSig = nullptr;
                    TGraphErrors* gSubPtSig = nullptr;
                    if (!subPtX.empty())
                    {
                        double yLoSig = 1e30, yHiSig = -1e30;
                        for (std::size_t is = 0; is < subPtSigY.size(); ++is)
                        {
                            yLoSig = std::min(yLoSig, subPtSigY[is] - subPtSigEY[is]);
                            yHiSig = std::max(yHiSig, subPtSigY[is] + subPtSigEY[is]);
                        }
                        const double sigPadPt = (yHiSig > yLoSig) ? 0.25 * (yHiSig - yLoSig) : 1.0;
                        
                        hFrSubPtSig = new TH1F(TString::Format("hFrSubSig_ptOv_%s_%s", H.variant.c_str(), cb.folder.c_str()).Data(),
                                               "", 100, subPtLoEdge, subPtHiEdge);
                        hFrSubPtSig->SetDirectory(nullptr);
                        hFrSubPtSig->SetStats(0);
                        hFrSubPtSig->SetMinimum(yLoSig - sigPadPt);
                        hFrSubPtSig->SetMaximum(yHiSig + sigPadPt);
                        hFrSubPtSig->GetYaxis()->SetTitle("#sigma^{Gauss}[GeV]");
                        hFrSubPtSig->GetYaxis()->SetTitleSize(0.14);
                        hFrSubPtSig->GetYaxis()->SetTitleOffset(0.32);
                        hFrSubPtSig->GetYaxis()->SetLabelSize(0.09);
                        hFrSubPtSig->GetYaxis()->SetLabelOffset(0.008);
                        hFrSubPtSig->GetYaxis()->SetNdivisions(505);
                        hFrSubPtSig->GetXaxis()->SetNdivisions(505);
                        hFrSubPtSig->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hFrSubPtSig->GetXaxis()->SetTitleSize(0.12);
                        hFrSubPtSig->GetXaxis()->SetTitleOffset(0.95);
                        hFrSubPtSig->GetXaxis()->SetLabelSize(0.10);
                        hFrSubPtSig->Draw();
                        
                        vector<double> exZeroSubPtSig(subPtX.size(), 0.0);
                        gSubPtSig = new TGraphErrors((int)subPtX.size(), &subPtX[0], &subPtSigY[0], &exZeroSubPtSig[0], &subPtSigEY[0]);
                        gSubPtSig->SetMarkerStyle(20);
                        gSubPtSig->SetMarkerSize(1.0);
                        gSubPtSig->SetMarkerColor(kBlack);
                        gSubPtSig->SetLineColor(kBlack);
                        gSubPtSig->SetLineWidth(2);
                        gSubPtSig->Draw("PE1 SAME");
                    }
                    
                    cPO.Modified();
                    cPO.Update();
                    SaveCanvas(cPO, JoinPath(centSubDir, "isoPtOverlay_counts.png"));
                    if (!perVariantIsoQAOnlyActive)
                        SaveCanvas(cPO, JoinPath(variantPtSummaryDir,
                                                 TString::Format("isoPtOverlay_counts_%s.png", cb.folder.c_str()).Data()));
                    
                    if (gSubPtSig) delete gSubPtSig;
                    if (hFrSubPtSig) delete hFrSubPtSig;
                    if (gSubPtMean) delete gSubPtMean;
                    if (hFrSubPtMean) delete hFrSubPtMean;
                    
                    for (auto* f : gaussCurvesP) delete f;
                    for (auto* h : hPts) delete h;
                }
            }
            
            // ── combined data + inclusive MC + photon+jet MC overlay with μ/σ subpanels ──
            if (!forEmbeddedSim && incMCvarTop && phoMCvarTop)
            {
                for (const auto& cb : centPlotBins)
                {
                    const string centDirComb = JoinPath(variantDir, cb.folder);
                    
                    // Pre-pass: accumulate Gaussian fit params across pT
                    vector<double> combDataX, combDataMuY, combDataMuEY, combDataSigY, combDataSigEY;
                    vector<double> combIncX,  combIncMuY,  combIncMuEY,  combIncSigY,  combIncSigEY;
                    vector<double> combPhoX,  combPhoMuY,  combPhoMuEY,  combPhoSigY,  combPhoSigEY;
                    
                    for (int iptPre = 0; iptPre < kNPtBins; ++iptPre)
                    {
                        const PtBin& bPre = PtBins()[iptPre];
                        if (bPre.lo < 10) continue;
                        const double ptC = 0.5 * (kPtEdges[(std::size_t)iptPre] + kPtEdges[(std::size_t)iptPre + 1]);
                        const string hName = "h_Eiso" + bPre.suffix + cb.suffix;
                        
                        auto FitPushComb = [&](TDirectory* dir, const string& hN,
                                               vector<double>& vx, vector<double>& vmy, vector<double>& vmey,
                                               vector<double>& vsy, vector<double>& vsey) {
                            TH1* hSrc = dynamic_cast<TH1*>(dir->Get(hN.c_str()));
                            if (!hSrc) return;
                            TH1* hTmp = CloneTH1(hSrc, TString::Format("hCombPre_%s_%d_%p", hN.c_str(), iptPre, (void*)dir).Data());
                            if (!hTmp) return;
                            hTmp->Rebin(10); EnsureSumw2(hTmp);
                            double gM, gS, gME, gSE;
                            if (FitGaussianIterative(hTmp, gM, gS, gME, gSE))
                            { vx.push_back(ptC); vmy.push_back(gM); vmey.push_back(gME); vsy.push_back(gS); vsey.push_back(gSE); }
                            delete hTmp;
                        };
                        
                        FitPushComb(aaTop,        hName, combDataX, combDataMuY, combDataMuEY, combDataSigY, combDataSigEY);
                        FitPushComb(incMCvarTop,  hName, combIncX,  combIncMuY,  combIncMuEY,  combIncSigY,  combIncSigEY);
                        FitPushComb(phoMCvarTop,  hName, combPhoX,  combPhoMuY,  combPhoMuEY,  combPhoSigY,  combPhoSigEY);
                    }
                    
                    // Compute μ range
                    double combMuLo = 1e30, combMuHi = -1e30;
                    auto CombUpdateMu = [&](const vector<double>& y, const vector<double>& ey) {
                        for (std::size_t i = 0; i < y.size(); ++i) {
                            combMuLo = std::min(combMuLo, y[i] - ey[i]);
                            combMuHi = std::max(combMuHi, y[i] + ey[i]); } };
                    CombUpdateMu(combDataMuY, combDataMuEY);
                    CombUpdateMu(combIncMuY,  combIncMuEY);
                    CombUpdateMu(combPhoMuY,  combPhoMuEY);
                    const bool combHaveSub = (combMuHi > combMuLo) &&
                    (!combDataX.empty() || !combIncX.empty() || !combPhoX.empty());
                    const double combMuPad = combHaveSub ? std::max(0.35 * (combMuHi - combMuLo), 0.5) : 1.0;
                    
                    // Compute σ range
                    double combSigLo = 1e30, combSigHi = -1e30;
                    auto CombUpdateSig = [&](const vector<double>& y, const vector<double>& ey) {
                        for (std::size_t i = 0; i < y.size(); ++i) {
                            combSigLo = std::min(combSigLo, y[i] - ey[i]);
                            combSigHi = std::max(combSigHi, y[i] + ey[i]); } };
                    CombUpdateSig(combDataSigY, combDataSigEY);
                    CombUpdateSig(combIncSigY,  combIncSigEY);
                    CombUpdateSig(combPhoSigY,  combPhoSigEY);
                    const double combSigPad = combHaveSub ? std::max(0.35 * (combSigHi - combSigLo), 0.5) : 1.0;
                    
                    for (const auto& b : PtBins())
                    {
                        if (b.lo < 10) continue;
                        const string ptDirComb = JoinPath(centDirComb, b.folder);
                        const string combOvDir = JoinPath(ptDirComb, "photonJetOverlays_inclusiveMCoverlays");
                        EnsureDir(combOvDir);
                        
                        const string hAAName = "h_Eiso" + b.suffix + cb.suffix;
                        TH1* hDataSrc = dynamic_cast<TH1*>(aaTop->Get(hAAName.c_str()));
                        TH1* hIncSrc  = dynamic_cast<TH1*>(incMCvarTop->Get(hAAName.c_str()));
                        TH1* hPhoSrc  = dynamic_cast<TH1*>(phoMCvarTop->Get(hAAName.c_str()));
                        if (!hDataSrc || !hIncSrc || !hPhoSrc) continue;
                        
                        TH1* hData = CloneTH1(hDataSrc, TString::Format("hData_combOv_%s_%s_%s", H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                        TH1* hInc  = CloneTH1(hIncSrc,  TString::Format("hInc_combOv_%s_%s_%s",  H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                        TH1* hPho  = CloneTH1(hPhoSrc,  TString::Format("hPho_combOv_%s_%s_%s",  H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                        if (!hData || !hInc || !hPho) { if (hData) delete hData; if (hInc) delete hInc; if (hPho) delete hPho; continue; }
                        
                        EnsureSumw2(hData); EnsureSumw2(hInc); EnsureSumw2(hPho);
                        hData->Rebin(10); hInc->Rebin(10); hPho->Rebin(10);
                        
                        const double intData = hData->Integral(0, hData->GetNbinsX() + 1);
                        const double intInc  = hInc->Integral(0, hInc->GetNbinsX() + 1);
                        const double intPho  = hPho->Integral(0, hPho->GetNbinsX() + 1);
                        if (!(intData > 0.0) || !(intInc > 0.0) || !(intPho > 0.0)) { delete hData; delete hInc; delete hPho; continue; }
                        
                        hData->Scale(1.0 / intData);
                        hInc->Scale(1.0 / intInc);
                        hPho->Scale(1.0 / intPho);
                        
                        // Data: black markers
                        hData->SetLineColor(kBlack);
                        hData->SetMarkerColor(kBlack);
                        hData->SetMarkerStyle(20);
                        hData->SetMarkerSize(1.0);
                        hData->SetLineWidth(2);
                        hData->SetFillStyle(0);
                        
                        // Photon+jet MC: red histogram
                        hPho->SetTitle("");
                        hPho->SetLineColor(kRed + 1);
                        hPho->SetLineWidth(2);
                        hPho->SetFillStyle(0);
                        hPho->SetMarkerSize(0.0);
                        TH1* hPhoForFit = (TH1*)hPho->Clone(TString::Format("hPhoFit_combOv_%s_%s_%s", H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                        hPhoForFit->SetDirectory(nullptr);
                        for (int ib = 0; ib <= hPho->GetNbinsX() + 1; ++ib) hPho->SetBinError(ib, 0.0);
                        
                        // Inclusive jet MC: blue histogram
                        hInc->SetLineColor(kBlue + 1);
                        hInc->SetLineWidth(2);
                        hInc->SetFillStyle(0);
                        hInc->SetMarkerSize(0.0);
                        TH1* hIncForFit = (TH1*)hInc->Clone(TString::Format("hIncFit_combOv_%s_%s_%s", H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                        hIncForFit->SetDirectory(nullptr);
                        for (int ib = 0; ib <= hInc->GetNbinsX() + 1; ++ib) hInc->SetBinError(ib, 0.0);
                        
                        const double ymx = std::max({hData->GetMaximum(), hInc->GetMaximum(), hPho->GetMaximum()});
                        
                        // Frame histogram for axis ownership (drawn first)
                        TH1* hFrame = (TH1*)hPho->Clone(TString::Format("hFr_combOv_%s_%s_%s", H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                        hFrame->SetDirectory(nullptr);
                        
                        TCanvas cCombOv(
                                        TString::Format("c_combOv_%s_%s_%s", H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                        "c_combOv", 900, combHaveSub ? 1000 : 700);
                        cCombOv.cd();
                        
                        // Upper pad
                        const double combPadLoEdge = combHaveSub ? 0.36 : 0.0;
                        TPad* padUpComb = new TPad("padUpComb", "padUpComb", 0.0, combPadLoEdge, 1.0, 1.0);
                        padUpComb->SetBottomMargin(combHaveSub ? 0.10 : 0.12);
                        padUpComb->SetLeftMargin(0.14);
                        padUpComb->SetRightMargin(0.04);
                        padUpComb->SetTopMargin(0.08);
                        padUpComb->Draw();
                        padUpComb->cd();
                        
                        hFrame->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                        hFrame->GetYaxis()->SetTitle("Normalized to unit area");
                        hFrame->GetXaxis()->SetTitleSize(0.045);
                        hFrame->GetYaxis()->SetTitleSize(0.055);
                        hFrame->GetXaxis()->SetLabelSize(0.045);
                        hFrame->GetYaxis()->SetLabelSize(0.045);
                        hFrame->GetYaxis()->SetTitleOffset(1.15);
                        hFrame->SetMinimum(0.0);
                        hFrame->SetMaximum((ymx > 0.0) ? (1.40 * ymx) : 1.0);
                        
                        hFrame->Draw("hist");
                        hInc->Draw("hist SAME");
                        hPho->Draw("hist SAME");
                        hData->Draw("E1 SAME");
                        
                        // Gaussian fit curves in matching colors
                        TF1* fDataGauss = DrawGaussFitCurve(hData, kBlack);
                        TF1* fPhoGauss  = DrawGaussFitCurve(hPhoForFit, kRed + 1);
                        TF1* fIncGauss  = DrawGaussFitCurve(hIncForFit, kBlue + 1);
                        
                        TLegend legComb(0.55, 0.65, 0.95, 0.85);
                        legComb.SetBorderSize(0);
                        legComb.SetFillStyle(0);
                        legComb.SetTextFont(42);
                        legComb.SetTextSize(0.038);
                        legComb.AddEntry(hData, TString::Format("AuAu data (%s)", H.label.c_str()).Data(), "ep");
                        legComb.AddEntry(hPho, "photon+jet embedded MC", "l");
                        legComb.AddEntry(hInc, "inclusive jet embedded MC", "l");
                        legComb.Draw();
                        
                        TLatex tCombTitle;
                        tCombTitle.SetNDC(true);
                        tCombTitle.SetTextFont(42);
                        tCombTitle.SetTextAlign(23);
                        tCombTitle.SetTextSize(0.038);
                        tCombTitle.DrawLatex(0.50, 0.97,
                                             "E_{T}^{iso} overlay: AuAu data vs photon+jet & inclusive embedded MC");
                        
                        TLatex tCombInfo;
                        tCombInfo.SetNDC(true);
                        tCombInfo.SetTextFont(42);
                        tCombInfo.SetTextAlign(13);
                        tCombInfo.SetTextSize(0.045);
                        tCombInfo.DrawLatex(0.22, 0.88, TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                        tCombInfo.DrawLatex(0.22, 0.82, TString::Format("p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                        
                        {
                            const bool combInfoRHS = (cb.lo >= 40);
                            const int   combInfoAlign = combInfoRHS ? 33 : 13;
                            const double combInfoX    = combInfoRHS ? 0.92 : 0.16;
                            const double combSphY1    = combInfoRHS ? 0.60 : 0.73;
                            TLatex tSph;
                            tSph.SetNDC(true);
                            tSph.SetTextFont(42);
                            tSph.SetTextAlign(combInfoAlign);
                            tSph.SetTextSize(0.042);
                            tSph.DrawLatex(combInfoX, combSphY1, "#bf{sPHENIX} #it{Internal}");
                            tSph.SetTextSize(0.034);
                            tSph.DrawLatex(combInfoX, combSphY1 - 0.06, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                        }
                        {
                            const bool combInfoRHS = (cb.lo >= 40);
                            const int   combInfoAlign = combInfoRHS ? 33 : 13;
                            const double combInfoX    = combInfoRHS ? 0.92 : 0.16;
                            const double combUEY1     = combInfoRHS ? 0.48 : 0.61;
                            TLatex tUE;
                            tUE.SetNDC(true);
                            tUE.SetTextFont(42);
                            tUE.SetTextAlign(combInfoAlign);
                            tUE.SetTextSize(0.030);
                            tUE.DrawLatex(combInfoX, combUEY1,        trigDisplayLabel.c_str());
                            tUE.DrawLatex(combInfoX, combUEY1 - 0.04, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                            tUE.DrawLatex(combInfoX, combUEY1 - 0.08, TString::Format("UE: %s", H.label.c_str()).Data());
                            tUE.DrawLatex(combInfoX, combUEY1 - 0.12, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        }
                        
                        // Middle pad: Gaussian mean vs pT
                        if (combHaveSub)
                        {
                            auto MakeSubGComb = [](const vector<double>& x, const vector<double>& y, const vector<double>& ey,
                                                   int marker, int color) -> TGraphErrors* {
                                if (x.empty()) return nullptr;
                                vector<double> ex(x.size(), 0.0);
                                TGraphErrors* g = new TGraphErrors((int)x.size(), &x[0], &y[0], &ex[0], &ey[0]);
                                g->SetMarkerStyle(marker); g->SetMarkerSize(1.0);
                                g->SetMarkerColor(color);  g->SetLineColor(color);
                                g->SetLineWidth(2); g->Draw("PE1 SAME"); return g; };
                            
                            cCombOv.cd();
                            TPad* padMidComb = new TPad("padMidComb", "padMidComb", 0.0, 0.21, 1.0, 0.36);
                            padMidComb->SetTopMargin(0.02);
                            padMidComb->SetBottomMargin(0.00);
                            padMidComb->SetLeftMargin(0.14);
                            padMidComb->SetRightMargin(0.04);
                            padMidComb->Draw();
                            padMidComb->cd();
                            
                            TH1F* hFrMuComb = new TH1F(
                                                       TString::Format("hFrMu_comb_%s_%s_%s", H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                                       "", 100, 10.0, kPtEdges.back());
                            hFrMuComb->SetDirectory(nullptr); hFrMuComb->SetStats(0);
                            hFrMuComb->SetMinimum(combMuLo - combMuPad);
                            hFrMuComb->SetMaximum(combMuHi + combMuPad);
                            hFrMuComb->GetYaxis()->SetTitle("#mu^{Gauss}[GeV]");
                            hFrMuComb->GetYaxis()->SetTitleSize(0.19);
                            hFrMuComb->GetYaxis()->SetTitleOffset(0.28);
                            hFrMuComb->GetYaxis()->SetLabelSize(0.14);
                            hFrMuComb->GetYaxis()->SetNdivisions(505);
                            hFrMuComb->GetXaxis()->SetTitle("");
                            hFrMuComb->GetXaxis()->SetTitleSize(0.0);
                            hFrMuComb->GetXaxis()->SetLabelSize(0.0);
                            hFrMuComb->GetXaxis()->SetTickLength(0.0);
                            hFrMuComb->Draw();
                            
                            TGraphErrors* gMuData = MakeSubGComb(combDataX, combDataMuY, combDataMuEY, 20, kBlack);
                            TGraphErrors* gMuPho  = MakeSubGComb(combPhoX,  combPhoMuY,  combPhoMuEY,  24, kRed+1);
                            TGraphErrors* gMuInc  = MakeSubGComb(combIncX,  combIncMuY,  combIncMuEY,  25, kBlue+1);
                            
                            // Bottom pad: Gaussian sigma vs pT
                            cCombOv.cd();
                            TPad* padBotComb = new TPad("padBotComb", "padBotComb", 0.0, 0.0, 1.0, 0.21);
                            padBotComb->SetTopMargin(0.00);
                            padBotComb->SetBottomMargin(0.30);
                            padBotComb->SetLeftMargin(0.14);
                            padBotComb->SetRightMargin(0.04);
                            padBotComb->Draw();
                            padBotComb->cd();
                            
                            TH1F* hFrSigComb = new TH1F(
                                                        TString::Format("hFrSig_comb_%s_%s_%s", H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                                        "", 100, 10.0, kPtEdges.back());
                            hFrSigComb->SetDirectory(nullptr); hFrSigComb->SetStats(0);
                            hFrSigComb->SetMinimum(combSigLo - combSigPad);
                            hFrSigComb->SetMaximum(combSigHi + combSigPad);
                            hFrSigComb->GetYaxis()->SetTitle("#sigma^{Gauss}[GeV]");
                            hFrSigComb->GetYaxis()->SetTitleSize(0.14);
                            hFrSigComb->GetYaxis()->SetTitleOffset(0.32);
                            hFrSigComb->GetYaxis()->SetLabelSize(0.10);
                            hFrSigComb->GetYaxis()->SetNdivisions(505);
                            hFrSigComb->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                            hFrSigComb->GetXaxis()->SetTitleSize(0.12);
                            hFrSigComb->GetXaxis()->SetTitleOffset(0.95);
                            hFrSigComb->GetXaxis()->SetLabelSize(0.10);
                            hFrSigComb->Draw();
                            
                            TGraphErrors* gSigData = MakeSubGComb(combDataX, combDataSigY, combDataSigEY, 20, kBlack);
                            TGraphErrors* gSigPho  = MakeSubGComb(combPhoX,  combPhoSigY,  combPhoSigEY,  24, kRed+1);
                            TGraphErrors* gSigInc  = MakeSubGComb(combIncX,  combIncSigY,  combIncSigEY,  25, kBlue+1);
                            
                            cCombOv.Modified();
                            cCombOv.Update();
                            
                            SaveCanvas(cCombOv, JoinPath(combOvDir, "Eiso_dataMC_overlay.png"));
                            
                            if (gSigData) delete gSigData;
                            if (gSigPho)  delete gSigPho;
                            if (gSigInc)  delete gSigInc;
                            delete hFrSigComb;
                            if (gMuData) delete gMuData;
                            if (gMuPho)  delete gMuPho;
                            if (gMuInc)  delete gMuInc;
                            delete hFrMuComb;
                        }
                        
                        if (!combHaveSub)
                            SaveCanvas(cCombOv, JoinPath(combOvDir, "Eiso_dataMC_overlay.png"));
                        
                        if (fDataGauss) delete fDataGauss;
                        if (fPhoGauss)  delete fPhoGauss;
                        if (fIncGauss)  delete fIncGauss;
                        delete hPhoForFit;
                        delete hIncForFit;
                        delete hFrame;
                        delete hData;
                        delete hInc;
                        delete hPho;
                    }
                    
                    // ── merged-pT combined overlay (photonJetOverlays_inclusiveMCoverlays) ──
                    {
                        struct MergedPtRange { int lo; int hi; };
                        const vector<MergedPtRange> mergedRanges = { {20, 35} };
                        
                        for (const auto& mr : mergedRanges)
                        {
                            vector<int> mBins;
                            for (int ip = 0; ip < kNPtBins; ++ip)
                            {
                                const PtBin& bp = PtBins()[ip];
                                if (bp.lo >= mr.lo && bp.hi <= mr.hi) mBins.push_back(ip);
                            }
                            if (mBins.empty()) continue;
                            
                            auto MergeCombHist = [&](TDirectory* dir, const string& centSuffix, const string& cloneName) -> TH1*
                            {
                                TH1* hSum = nullptr;
                                for (int ip : mBins)
                                {
                                    const PtBin& bp = PtBins()[ip];
                                    const string hName = "h_Eiso" + bp.suffix + centSuffix;
                                    TH1* hSrc = dynamic_cast<TH1*>(dir->Get(hName.c_str()));
                                    if (!hSrc) continue;
                                    if (!hSum) { hSum = CloneTH1(hSrc, cloneName.c_str()); if (hSum) EnsureSumw2(hSum); }
                                    else       { hSum->Add(hSrc); }
                                }
                                return hSum;
                            };
                            
                            const string mFolder = TString::Format("pT_%d_%d", mr.lo, mr.hi).Data();
                            const string mPtDir  = JoinPath(centDirComb, mFolder);
                            const string mCombDir = JoinPath(mPtDir, "photonJetOverlays_inclusiveMCoverlays");
                            EnsureDir(mCombDir);
                            
                            const string mTag = H.variant + "_" + cb.folder + "_" + mFolder;
                            TH1* hDataM = MergeCombHist(aaTop,       cb.suffix, "hData_combMerged_" + mTag);
                            TH1* hIncM  = MergeCombHist(incMCvarTop, cb.suffix, "hInc_combMerged_"  + mTag);
                            TH1* hPhoM  = MergeCombHist(phoMCvarTop, cb.suffix, "hPho_combMerged_"  + mTag);
                            if (!hDataM || !hIncM || !hPhoM) { if (hDataM) delete hDataM; if (hIncM) delete hIncM; if (hPhoM) delete hPhoM; continue; }
                            
                            EnsureSumw2(hDataM); EnsureSumw2(hIncM); EnsureSumw2(hPhoM);
                            hDataM->Rebin(10); hIncM->Rebin(10); hPhoM->Rebin(10);
                            
                            const double intD = hDataM->Integral(0, hDataM->GetNbinsX()+1);
                            const double intI = hIncM->Integral(0,  hIncM->GetNbinsX()+1);
                            const double intP = hPhoM->Integral(0,  hPhoM->GetNbinsX()+1);
                            if (!(intD > 0.0) || !(intI > 0.0) || !(intP > 0.0)) { delete hDataM; delete hIncM; delete hPhoM; continue; }
                            
                            hDataM->Scale(1.0/intD); hIncM->Scale(1.0/intI); hPhoM->Scale(1.0/intP);
                            
                            hDataM->SetLineColor(kBlack); hDataM->SetMarkerColor(kBlack); hDataM->SetMarkerStyle(20);
                            hDataM->SetMarkerSize(1.0); hDataM->SetLineWidth(2); hDataM->SetFillStyle(0);
                            hPhoM->SetTitle(""); hPhoM->SetLineColor(kRed+1); hPhoM->SetLineWidth(2);
                            hPhoM->SetFillStyle(0); hPhoM->SetMarkerSize(0.0);
                            TH1* hPhoMForFit = (TH1*)hPhoM->Clone(("hPhoMFit_" + mTag).c_str());
                            hPhoMForFit->SetDirectory(nullptr);
                            for (int ib=0; ib<=hPhoM->GetNbinsX()+1; ++ib) hPhoM->SetBinError(ib,0.0);
                            hIncM->SetLineColor(kBlue+1); hIncM->SetLineWidth(2);
                            hIncM->SetFillStyle(0); hIncM->SetMarkerSize(0.0);
                            TH1* hIncMForFit = (TH1*)hIncM->Clone(("hIncMFit_" + mTag).c_str());
                            hIncMForFit->SetDirectory(nullptr);
                            for (int ib=0; ib<=hIncM->GetNbinsX()+1; ++ib) hIncM->SetBinError(ib,0.0);
                            
                            const double ymxM = std::max({hDataM->GetMaximum(), hIncM->GetMaximum(), hPhoM->GetMaximum()});
                            TH1* hFrM = (TH1*)hPhoM->Clone(("hFr_combMerged_" + mTag).c_str());
                            hFrM->SetDirectory(nullptr);
                            
                            const int mCanvasH = combHaveSub ? 1000 : 700;
                            TCanvas cM(("c_combMerged_" + mTag).c_str(), "c_combMerged", 900, mCanvasH);
                            cM.cd();
                            
                            const double mPadLoEdge = combHaveSub ? 0.36 : 0.0;
                            TPad* padUpM = new TPad("padUpM", "padUpM", 0.0, mPadLoEdge, 1.0, 1.0);
                            padUpM->SetBottomMargin(combHaveSub ? 0.10 : 0.12);
                            padUpM->SetLeftMargin(0.14); padUpM->SetRightMargin(0.04); padUpM->SetTopMargin(0.08);
                            padUpM->Draw(); padUpM->cd();
                            
                            hFrM->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                            hFrM->GetYaxis()->SetTitle("Normalized to unit area");
                            hFrM->GetXaxis()->SetTitleSize(0.045); hFrM->GetYaxis()->SetTitleSize(0.055);
                            hFrM->GetXaxis()->SetLabelSize(0.045); hFrM->GetYaxis()->SetLabelSize(0.045);
                            hFrM->GetYaxis()->SetTitleOffset(1.15);
                            hFrM->SetMinimum(0.0); hFrM->SetMaximum((ymxM > 0.0) ? (1.40*ymxM) : 1.0);
                            hFrM->Draw("hist"); hIncM->Draw("hist SAME"); hPhoM->Draw("hist SAME"); hDataM->Draw("E1 SAME");
                            
                            TF1* fDG = DrawGaussFitCurve(hDataM, kBlack);
                            TF1* fPG = DrawGaussFitCurve(hPhoMForFit, kRed+1);
                            TF1* fIG = DrawGaussFitCurve(hIncMForFit, kBlue+1);
                            
                            TLegend lgM(0.55, 0.68, 0.95, 0.88); lgM.SetBorderSize(0); lgM.SetFillStyle(0); lgM.SetTextFont(42); lgM.SetTextSize(0.034);
                            lgM.AddEntry(hDataM, TString::Format("AuAu data (%s)", H.label.c_str()).Data(), "ep");
                            lgM.AddEntry(hPhoM,  "photon+jet embedded MC", "l");
                            lgM.AddEntry(hIncM,  "inclusive jet embedded MC", "l");
                            lgM.Draw();
                            
                            TLatex tMt; tMt.SetNDC(true); tMt.SetTextFont(42); tMt.SetTextAlign(23); tMt.SetTextSize(0.038);
                            tMt.DrawLatex(0.50, 0.97, "E_{T}^{iso} overlay: AuAu data vs photon+jet & inclusive embedded MC");
                            TLatex tMi; tMi.SetNDC(true); tMi.SetTextFont(42); tMi.SetTextAlign(13); tMi.SetTextSize(0.045);
                            tMi.DrawLatex(0.22, 0.88, TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                            tMi.DrawLatex(0.22, 0.82, TString::Format("p_{T}^{#gamma} = %d-%d GeV", mr.lo, mr.hi).Data());
                            {
                                const bool mInfoRHS = (cb.lo >= 20);
                                const int   mInfoAlign = mInfoRHS ? 33 : 13;
                                const double mInfoX    = mInfoRHS ? 0.92 : 0.16;
                                const double mSphY1    = mInfoRHS ? 0.60 : 0.73;
                                TLatex tS; tS.SetNDC(true); tS.SetTextFont(42); tS.SetTextAlign(mInfoAlign);
                                tS.SetTextSize(0.042);
                                tS.DrawLatex(mInfoX, mSphY1, "#bf{sPHENIX} #it{Internal}");
                                tS.SetTextSize(0.034);
                                tS.DrawLatex(mInfoX, mSphY1 - 0.06, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                            }
                            {
                                const bool mInfoRHS = (cb.lo >= 20);
                                const int   mInfoAlign = mInfoRHS ? 33 : 13;
                                const double mInfoX    = mInfoRHS ? 0.92 : 0.16;
                                const double mUEY1     = mInfoRHS ? 0.48 : 0.61;
                                TLatex tU; tU.SetNDC(true); tU.SetTextFont(42); tU.SetTextAlign(mInfoAlign);
                                tU.SetTextSize(0.026);
                                tU.DrawLatex(mInfoX, mUEY1,        trigDisplayLabel.c_str());
                                tU.SetTextSize(0.030);
                                tU.DrawLatex(mInfoX, mUEY1 - 0.04, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                                tU.DrawLatex(mInfoX, mUEY1 - 0.08, TString::Format("UE: %s", H.label.c_str()).Data());
                                tU.DrawLatex(mInfoX, mUEY1 - 0.12, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                            }
                            
                            // Middle pad: Gaussian mean vs pT
                            if (combHaveSub)
                            {
                                auto MakeSubGM = [](const vector<double>& x, const vector<double>& y, const vector<double>& ey,
                                                    int marker, int color) -> TGraphErrors* {
                                    if (x.empty()) return nullptr;
                                    vector<double> ex(x.size(), 0.0);
                                    TGraphErrors* g = new TGraphErrors((int)x.size(), &x[0], &y[0], &ex[0], &ey[0]);
                                    g->SetMarkerStyle(marker); g->SetMarkerSize(1.0);
                                    g->SetMarkerColor(color); g->SetLineColor(color);
                                    g->SetLineWidth(2); g->Draw("PE1 SAME"); return g; };
                                
                                cM.cd();
                                TPad* padMidM = new TPad("padMidM", "padMidM", 0.0, 0.21, 1.0, 0.36);
                                padMidM->SetTopMargin(0.02); padMidM->SetBottomMargin(0.00);
                                padMidM->SetLeftMargin(0.14); padMidM->SetRightMargin(0.04);
                                padMidM->Draw(); padMidM->cd();
                                
                                TH1F* hFrMuM = new TH1F(("hFrMu_combM_" + mTag).c_str(), "", 100, 10.0, kPtEdges.back());
                                hFrMuM->SetDirectory(nullptr); hFrMuM->SetStats(0);
                                hFrMuM->SetMinimum(combMuLo - combMuPad); hFrMuM->SetMaximum(combMuHi + combMuPad);
                                hFrMuM->GetYaxis()->SetTitle("#mu^{Gauss}[GeV]");
                                hFrMuM->GetYaxis()->SetTitleSize(0.19); hFrMuM->GetYaxis()->SetTitleOffset(0.28);
                                hFrMuM->GetYaxis()->SetLabelSize(0.14); hFrMuM->GetYaxis()->SetNdivisions(505);
                                hFrMuM->GetXaxis()->SetTitle(""); hFrMuM->GetXaxis()->SetTitleSize(0.0);
                                hFrMuM->GetXaxis()->SetLabelSize(0.0); hFrMuM->GetXaxis()->SetTickLength(0.0);
                                hFrMuM->Draw();
                                
                                TGraphErrors* gMuDataM = MakeSubGM(combDataX, combDataMuY, combDataMuEY, 20, kBlack);
                                TGraphErrors* gMuPhoM  = MakeSubGM(combPhoX,  combPhoMuY,  combPhoMuEY,  24, kRed+1);
                                TGraphErrors* gMuIncM  = MakeSubGM(combIncX,  combIncMuY,  combIncMuEY,  25, kBlue+1);
                                
                                // Bottom pad: Gaussian sigma vs pT
                                cM.cd();
                                TPad* padBotM = new TPad("padBotM", "padBotM", 0.0, 0.0, 1.0, 0.21);
                                padBotM->SetTopMargin(0.00); padBotM->SetBottomMargin(0.30);
                                padBotM->SetLeftMargin(0.14); padBotM->SetRightMargin(0.04);
                                padBotM->Draw(); padBotM->cd();
                                
                                TH1F* hFrSigM = new TH1F(("hFrSig_combM_" + mTag).c_str(), "", 100, 10.0, kPtEdges.back());
                                hFrSigM->SetDirectory(nullptr); hFrSigM->SetStats(0);
                                hFrSigM->SetMinimum(combSigLo - combSigPad); hFrSigM->SetMaximum(combSigHi + combSigPad);
                                hFrSigM->GetYaxis()->SetTitle("#sigma^{Gauss}[GeV]");
                                hFrSigM->GetYaxis()->SetTitleSize(0.14); hFrSigM->GetYaxis()->SetTitleOffset(0.32);
                                hFrSigM->GetYaxis()->SetLabelSize(0.10); hFrSigM->GetYaxis()->SetNdivisions(505);
                                hFrSigM->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                                hFrSigM->GetXaxis()->SetTitleSize(0.12); hFrSigM->GetXaxis()->SetTitleOffset(0.95);
                                hFrSigM->GetXaxis()->SetLabelSize(0.10);
                                hFrSigM->Draw();
                                
                                TGraphErrors* gSigDataM = MakeSubGM(combDataX, combDataSigY, combDataSigEY, 20, kBlack);
                                TGraphErrors* gSigPhoM  = MakeSubGM(combPhoX,  combPhoSigY,  combPhoSigEY,  24, kRed+1);
                                TGraphErrors* gSigIncM  = MakeSubGM(combIncX,  combIncSigY,  combIncSigEY,  25, kBlue+1);
                                
                                cM.Modified(); cM.Update();
                                SaveCanvas(cM, JoinPath(mCombDir, "Eiso_dataMC_overlay.png"));
                                
                                if (gSigDataM) delete gSigDataM; if (gSigPhoM) delete gSigPhoM; if (gSigIncM) delete gSigIncM;
                                delete hFrSigM;
                                if (gMuDataM) delete gMuDataM; if (gMuPhoM) delete gMuPhoM; if (gMuIncM) delete gMuIncM;
                                delete hFrMuM;
                            }
                            
                            if (!combHaveSub)
                                SaveCanvas(cM, JoinPath(mCombDir, "Eiso_dataMC_overlay.png"));
                            
                            if (fDG) delete fDG; if (fPG) delete fPG; if (fIG) delete fIG;
                            delete hPhoMForFit; delete hIncMForFit;
                            delete hFrM; delete hDataM; delete hIncM; delete hPhoM;
                        }
                    }
                }
            }
            
            if (!forEmbeddedSim)
            {
                const int centPtColors[] = {kBlack, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 1,
                                            kCyan + 1, kYellow + 2, kViolet + 1};
                
                // --- direct output to variantDir: Gaussian mean vs pT, centrality overlay + pp ---
                {
                    vector<TGraphErrors*> gCentMean;
                    vector<string>        centMeanLabels;
                    vector<double>        xPPMean, exPPMean, yPPMean, eyPPMean;
                    double yMinMean = std::numeric_limits<double>::max();
                    double yMaxMean = -std::numeric_limits<double>::max();
                    
                    for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                    {
                        const auto& cb = centBins[ic];
                        vector<double> xP, exP, yP, eyP;
                        
                        for (int ipt = 0; ipt < kNPtBins; ++ipt)
                        {
                            if (!gaussFilled[ivH][ipt][ic]) continue;
                            
                            xP.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                            exP.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                            yP.push_back(gaussMean[ivH][ipt][ic]);
                            eyP.push_back(gaussMeanErr[ivH][ipt][ic]);
                            
                            yMinMean = std::min(yMinMean, yP.back() - eyP.back());
                            yMaxMean = std::max(yMaxMean, yP.back() + eyP.back());
                        }
                        
                        if (xP.empty()) continue;
                        
                        TGraphErrors* g = new TGraphErrors((int)xP.size(), &xP[0], &yP[0], nullptr, &eyP[0]);
                        const int ci = (ic < 8) ? centPtColors[ic] : kBlack;
                        g->SetLineWidth(2);
                        g->SetLineColor(ci);
                        g->SetMarkerColor(ci);
                        g->SetMarkerStyle(20);
                        g->SetMarkerSize(1.2);
                        gCentMean.push_back(g);
                        centMeanLabels.push_back(TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                    }
                    
                    for (int ipt = 0; ipt < kNPtBins; ++ipt)
                    {
                        if (!ppGaussFilledByPt[ipt]) continue;
                        
                        xPPMean.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                        exPPMean.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                        yPPMean.push_back(ppGaussMeanByPt[ipt]);
                        eyPPMean.push_back(ppGaussMeanErrByPt[ipt]);
                        
                        yMinMean = std::min(yMinMean, yPPMean.back() - eyPPMean.back());
                        yMaxMean = std::max(yMaxMean, yPPMean.back() + eyPPMean.back());
                    }
                    
                    if (!gCentMean.empty() || !xPPMean.empty())
                    {
                        if (!std::isfinite(yMinMean) || !std::isfinite(yMaxMean))
                        {
                            yMinMean = 0.0;
                            yMaxMean = 1.0;
                        }
                        const double padMean = (yMaxMean > yMinMean) ? (0.20 * (yMaxMean - yMinMean)) : 0.5;
                        
                        TCanvas cMeanPtCent(
                            TString::Format("c_gaussMean_ppCentOverlay_vsPt_%s_%s",
                                            trigAA.c_str(), H.variant.c_str()).Data(),
                            "c_gaussMean_ppCentOverlay_vsPt", 900, 700);
                        ApplyCanvasMargins1D(cMeanPtCent);
                        cMeanPtCent.cd();
                        
                        TH1F hFrameMean(
                            TString::Format("hFrame_gaussMean_ppCentOverlay_vsPt_%s_%s",
                                            trigAA.c_str(), H.variant.c_str()).Data(),
                            "", 100, kPtEdges.front(), kPtEdges.back());
                        hFrameMean.SetDirectory(nullptr);
                        hFrameMean.SetStats(0);
                        hFrameMean.SetMinimum(yMinMean - padMean);
                        hFrameMean.SetMaximum(yMaxMean + padMean);
                        hFrameMean.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hFrameMean.GetYaxis()->SetTitle("#mu^{Gauss} [GeV]");
                        hFrameMean.GetXaxis()->SetTitleSize(0.055);
                        hFrameMean.GetYaxis()->SetTitleSize(0.055);
                        hFrameMean.GetXaxis()->SetLabelSize(0.045);
                        hFrameMean.GetYaxis()->SetLabelSize(0.045);
                        hFrameMean.GetYaxis()->SetTitleOffset(1.15);
                        hFrameMean.Draw();
                        
                        TGraphErrors* gPPMean = nullptr;
                        if (!xPPMean.empty())
                        {
                            gPPMean = new TGraphErrors((int)xPPMean.size(), &xPPMean[0], &yPPMean[0], nullptr, &eyPPMean[0]);
                            gPPMean->SetLineWidth(2);
                            gPPMean->SetLineColor(kRed + 1);
                            gPPMean->SetMarkerColor(kRed + 1);
                            gPPMean->SetMarkerStyle(24);
                            gPPMean->SetMarkerSize(1.3);
                            gPPMean->Draw("PE1 SAME");
                        }
                        
                        for (auto* g : gCentMean) g->Draw("PE1 SAME");
                        
                        TLegend legMean(0.19, 0.74, 0.47, 0.90);
                        legMean.SetBorderSize(0);
                        legMean.SetFillStyle(0);
                        legMean.SetTextFont(42);
                        legMean.SetTextSize(0.026);
                        legMean.SetNColumns(2);
                        if (gPPMean) legMean.AddEntry(gPPMean, "pp", "ep");
                        for (std::size_t ig = 0; ig < gCentMean.size(); ++ig)
                            legMean.AddEntry(gCentMean[ig], centMeanLabels[ig].c_str(), "ep");
                        legMean.Draw();
                        
                        TLatex tMeanTitle;
                        tMeanTitle.SetNDC(true);
                        tMeanTitle.SetTextFont(42);
                        tMeanTitle.SetTextAlign(23);
                        tMeanTitle.SetTextSize(0.042);
                        tMeanTitle.DrawLatex(0.50, 0.98,
                                             TString::Format("Gaussian #mu vs p_{T}^{#gamma}, centrality overlay + pp (%s)", H.label.c_str()).Data());
                        
                        TLatex tMeanInfo;
                        tMeanInfo.SetNDC(true);
                        tMeanInfo.SetTextFont(42);
                        tMeanInfo.SetTextAlign(33);
                        tMeanInfo.SetTextSize(0.028);
                        tMeanInfo.DrawLatex(0.94, 0.88, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
                        tMeanInfo.DrawLatex(0.94, 0.84, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                        tMeanInfo.DrawLatex(0.94, 0.80, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        tMeanInfo.DrawLatex(0.94, 0.76, TString::Format("UE: %s", H.label.c_str()).Data());
                        
                        SaveCanvas(cMeanPtCent, JoinPath(variantDir, "gaussMean_pp_auauCentralityOverlay_vs_pT.png"));
                        
                        if (gPPMean) delete gPPMean;
                    }
                    
                    for (auto* g : gCentMean) delete g;
                }
                
                // --- direct output to variantDir: Gaussian sigma vs pT, centrality overlay + pp ---
                {
                    vector<TGraphErrors*> gCentSigma;
                    vector<string>        centSigmaLabels;
                    vector<double>        xPPSigma, exPPSigma, yPPSigma, eyPPSigma;
                    double yMinSigma = std::numeric_limits<double>::max();
                    double yMaxSigma = -std::numeric_limits<double>::max();
                    
                    for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                    {
                        const auto& cb = centBins[ic];
                        vector<double> xP, exP, yP, eyP;
                        
                        for (int ipt = 0; ipt < kNPtBins; ++ipt)
                        {
                            if (!gaussFilled[ivH][ipt][ic]) continue;
                            
                            xP.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                            exP.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                            yP.push_back(gaussSigma[ivH][ipt][ic]);
                            eyP.push_back(gaussSigmaErr[ivH][ipt][ic]);
                            
                            yMinSigma = std::min(yMinSigma, yP.back() - eyP.back());
                            yMaxSigma = std::max(yMaxSigma, yP.back() + eyP.back());
                        }
                        
                        if (xP.empty()) continue;
                        
                        
                        TGraphErrors* g = new TGraphErrors((int)xP.size(), &xP[0], &yP[0], nullptr, &eyP[0]);
                        const int ci = (ic < 8) ? centPtColors[ic] : kBlack;
                        g->SetLineWidth(2);
                        g->SetLineColor(ci);
                        g->SetMarkerColor(ci);
                        g->SetMarkerStyle(20);
                        g->SetMarkerSize(1.2);
                        gCentSigma.push_back(g);
                        centSigmaLabels.push_back(TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                    }
                    
                    for (int ipt = 0; ipt < kNPtBins; ++ipt)
                    {
                        if (!ppGaussFilledByPt[ipt]) continue;
                        
                        xPPSigma.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                        exPPSigma.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                        yPPSigma.push_back(ppGaussSigmaByPt[ipt]);
                        eyPPSigma.push_back(ppGaussSigmaErrByPt[ipt]);
                        
                        yMinSigma = std::min(yMinSigma, yPPSigma.back() - eyPPSigma.back());
                        yMaxSigma = std::max(yMaxSigma, yPPSigma.back() + eyPPSigma.back());
                    }
                    
                    if (!gCentSigma.empty() || !xPPSigma.empty())
                    {
                        if (!std::isfinite(yMinSigma) || !std::isfinite(yMaxSigma))
                        {
                            yMinSigma = 0.0;
                            yMaxSigma = 1.0;
                        }
                        const double padSigma = (yMaxSigma > yMinSigma) ? (0.20 * (yMaxSigma - yMinSigma)) : 0.5;
                        
                        TCanvas cSigmaPtCent(
                            TString::Format("c_gaussSigma_ppCentOverlay_vsPt_%s_%s",
                                            trigAA.c_str(), H.variant.c_str()).Data(),
                            "c_gaussSigma_ppCentOverlay_vsPt", 900, 700);
                        ApplyCanvasMargins1D(cSigmaPtCent);
                        cSigmaPtCent.cd();
                        
                        TH1F hFrameSigma(
                            TString::Format("hFrame_gaussSigma_ppCentOverlay_vsPt_%s_%s",
                                            trigAA.c_str(), H.variant.c_str()).Data(),
                            "", 100, kPtEdges.front(), kPtEdges.back());
                        hFrameSigma.SetDirectory(nullptr);
                        hFrameSigma.SetStats(0);
                        hFrameSigma.SetMinimum(0.0);
                        hFrameSigma.SetMaximum(1.60 * std::max(yMaxSigma, 0.1));
                        hFrameSigma.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hFrameSigma.GetYaxis()->SetTitle("#sigma^{Gauss} [GeV]");
                        hFrameSigma.GetXaxis()->SetTitleSize(0.055);
                        hFrameSigma.GetYaxis()->SetTitleSize(0.055);
                        hFrameSigma.GetXaxis()->SetLabelSize(0.045);
                        hFrameSigma.GetYaxis()->SetLabelSize(0.045);
                        hFrameSigma.GetYaxis()->SetTitleOffset(1.15);
                        hFrameSigma.Draw();
                        
                        TGraphErrors* gPPSigma = nullptr;
                        if (!xPPSigma.empty())
                        {
                            gPPSigma = new TGraphErrors((int)xPPSigma.size(), &xPPSigma[0], &yPPSigma[0], nullptr, &eyPPSigma[0]);
                            gPPSigma->SetLineWidth(2);
                            gPPSigma->SetLineColor(kRed + 1);
                            gPPSigma->SetMarkerColor(kRed + 1);
                            gPPSigma->SetMarkerStyle(24);
                            gPPSigma->SetMarkerSize(1.3);
                            gPPSigma->Draw("PE1 SAME");
                        }
                        
                        for (auto* g : gCentSigma) g->Draw("PE1 SAME");
                        
                        TLegend legSigma(0.19, 0.74, 0.47, 0.90);
                        legSigma.SetBorderSize(0);
                        legSigma.SetFillStyle(0);
                        legSigma.SetTextFont(42);
                        legSigma.SetTextSize(0.026);
                        legSigma.SetNColumns(2);
                        if (gPPSigma) legSigma.AddEntry(gPPSigma, "pp", "ep");
                        for (std::size_t ig = 0; ig < gCentSigma.size(); ++ig)
                            legSigma.AddEntry(gCentSigma[ig], centSigmaLabels[ig].c_str(), "ep");
                        legSigma.Draw();
                        
                        TLatex tSigmaTitle;
                        tSigmaTitle.SetNDC(true);
                        tSigmaTitle.SetTextFont(42);
                        tSigmaTitle.SetTextAlign(23);
                        tSigmaTitle.SetTextSize(0.042);
                        tSigmaTitle.DrawLatex(0.50, 0.98,
                                              TString::Format("Gaussian #sigma vs p_{T}^{#gamma}, centrality overlay + pp (%s)", H.label.c_str()).Data());
                        
                        TLatex tSigmaInfo;
                        tSigmaInfo.SetNDC(true);
                        tSigmaInfo.SetTextFont(42);
                        tSigmaInfo.SetTextAlign(33);
                        tSigmaInfo.SetTextSize(0.028);
                        tSigmaInfo.DrawLatex(0.94, 0.88, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
                        tSigmaInfo.DrawLatex(0.94, 0.84, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                        tSigmaInfo.DrawLatex(0.94, 0.80, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        tSigmaInfo.DrawLatex(0.94, 0.76, TString::Format("UE: %s", H.label.c_str()).Data());
                        
                        SaveCanvas(cSigmaPtCent, JoinPath(variantDir, "gaussSigma_pp_auauCentralityOverlay_vs_pT.png"));
                        
                        if (gPPSigma) delete gPPSigma;
                    }
                    
                    for (auto* g : gCentSigma) delete g;
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
        
        if (!generateISOpTcentOverlaysONLY && !perVariantIsoQAOnlyActive)
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
        } // end if (!generateISOpTcentOverlaysONLY && !perVariantIsoQAOnlyActive) — perVariantOverlay + multi-variant sections
        
        if (!forEmbeddedSim)
        {
            
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
                            const double padS = (yMaxS > yMinS) ? 0.30 * (yMaxS - yMinS) : 0.5;
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
                            TLegend lgGS(0.78, 0.68, 0.86, 0.92);
                            lgGS.SetBorderSize(0); lgGS.SetFillStyle(0); lgGS.SetTextFont(42); lgGS.SetTextSize(0.032);
                            for (std::size_t ig = 0; ig < sGraphs.size(); ++ig)
                                lgGS.AddEntry(sGraphs[ig], handles[sIndices[ig]].label.c_str(), "ep");
                            lgGS.Draw();
                            TLatex tGS; tGS.SetNDC(true); tGS.SetTextFont(42); tGS.SetTextAlign(23); tGS.SetTextSize(0.042);
                            tGS.DrawLatex(0.50, 0.98,
                                          TString::Format("Gaussian #sigma vs Centrality, p_{T}^{#gamma} %d-%d GeV", b.lo, b.hi).Data());
                            TLatex tGSi; tGSi.SetNDC(true); tGSi.SetTextFont(42); tGSi.SetTextAlign(13); tGSi.SetTextSize(0.028);
                            tGSi.DrawLatex(0.20, 0.89, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
                            tGSi.DrawLatex(0.20, 0.85, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
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
                
                // -- baseVariant only: PPG12-style isolation cutoff fits from photon+jet embedded SIM --
                if (!forEmbeddedSim)
                {
                    const string baseVariantPtSummaryDir = JoinPath(ptSummaryBase, "baseVariant");
                    EnsureDir(baseVariantPtSummaryDir);
                    
                    const string phoCutIn = ResolvePhotonJetVariantInput("baseVariant");
                    if (!phoCutIn.empty())
                    {
                        TFile* fPhoCut = TFile::Open(phoCutIn.c_str(), "READ");
                        if (fPhoCut && !fPhoCut->IsZombie())
                        {
                            TDirectory* phoCutTop = fPhoCut->GetDirectory(kDirSIM.c_str());
                            if (!phoCutTop) phoCutTop = fPhoCut;
                            
                            auto FindEfficiencyCut = [&](TH1* hIn, double eff, double& cut, double& cutErr) -> bool
                            {
                                cut = 0.0;
                                cutErr = 0.0;
                                if (!hIn) return false;
                                
                                const int nb = hIn->GetNbinsX();
                                const double total = hIn->Integral(1, nb);
                                if (!(total > 0.0)) return false;
                                
                                double running = 0.0;
                                int ibCut = nb;
                                for (int ib = 1; ib <= nb; ++ib)
                                {
                                    running += hIn->GetBinContent(ib);
                                    if ((running / total) >= eff)
                                    {
                                        ibCut = ib;
                                        break;
                                    }
                                }
                                
                                const double binLo = hIn->GetXaxis()->GetBinLowEdge(ibCut);
                                const double binHi = hIn->GetXaxis()->GetBinUpEdge(ibCut);
                                const double prev = running - hIn->GetBinContent(ibCut);
                                const double binC = hIn->GetBinContent(ibCut);
                                const double target = eff * total;
                                
                                if (binC > 0.0)
                                {
                                    const double frac = std::min(1.0, std::max(0.0, (target - prev) / binC));
                                    cut = binLo + frac * (binHi - binLo);
                                }
                                else
                                {
                                    cut = 0.5 * (binLo + binHi);
                                }
                                
                                cutErr = 0.5 * (binHi - binLo);
                                return std::isfinite(cut);
                            };
                            
                            const string embeddedLabel = PhotonEmbeddedShortLabel().empty()
                            ? "Photon+Jet Embedded SIM"
                            : ("Photon+Jet " + PhotonEmbeddedShortLabel() + " Embedded SIM");
                            
                            // Open matching pp SIM (non-embedded) for PPG12 overlay
                            TFile* fMatchedPPSim = nullptr;
                            TDirectory* matchedPPSimTop = nullptr;
                            string matchedPPLabel;
                            {
                                SimSample matchedSample = SimSample::kPhotonJet20; // default
                                if (bothPhoton10and20simEmbedded)
                                {
                                    matchedSample = SimSample::kPhotonJet5And10And20Merged;
                                    matchedPPLabel = "Pythia 5+10+20";
                                }
                                else if (isPhotonJet10Embedded)
                                {
                                    matchedSample = SimSample::kPhotonJet10;
                                    matchedPPLabel = "Pythia 10";
                                }
                                else if (isPhotonJet20Embedded)
                                {
                                    matchedSample = SimSample::kPhotonJet20;
                                    matchedPPLabel = "Pythia 20";
                                }
                                else
                                {
                                    matchedPPLabel = "Pythia 20";
                                }
                                const string matchedPath = SimInputPathForSample(matchedSample);
                                if (!matchedPath.empty())
                                {
                                    fMatchedPPSim = TFile::Open(matchedPath.c_str(), "READ");
                                    if (fMatchedPPSim && !fMatchedPPSim->IsZombie())
                                        matchedPPSimTop = fMatchedPPSim->GetDirectory(kDirSIM.c_str());
                                    else { if (fMatchedPPSim) { fMatchedPPSim->Close(); delete fMatchedPPSim; fMatchedPPSim = nullptr; } }
                                }
                            }
                            
                            // Derive embedded legend tag from active sample
                            const string embLegTag = PhotonEmbeddedShortLabel().empty()
                            ? "Pythia Emb" : ("Pythia Emb " + PhotonEmbeddedShortLabel());
                            
                            for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                            {
                                const auto& cb = centBins[ic];
                                
                                vector<double> xCut;
                                vector<double> exCut;
                                vector<double> y70;
                                vector<double> ey70;
                                vector<double> y80;
                                vector<double> ey80;
                                vector<double> y90;
                                vector<double> ey90;
                                
                                double yMinEff = std::numeric_limits<double>::max();
                                double yMaxEff = -std::numeric_limits<double>::max();
                                
                                for (int ipt = 1; ipt < kNPtBins; ++ipt)
                                {
                                    const PtBin& b = PtBins()[ipt];
                                    const string hSigName = "h_EisoReco_truthSigMatched" + b.suffix + cb.suffix;
                                    TH1* hSigSrc = dynamic_cast<TH1*>(phoCutTop->Get(hSigName.c_str()));
                                    if (!hSigSrc || hSigSrc->GetEntries() <= 0.0) continue;
                                    
                                    double cut70 = 0.0, err70 = 0.0;
                                    double cut80 = 0.0, err80 = 0.0;
                                    double cut90 = 0.0, err90 = 0.0;
                                    
                                    if (!FindEfficiencyCut(hSigSrc, 0.70, cut70, err70)) continue;
                                    if (!FindEfficiencyCut(hSigSrc, 0.80, cut80, err80)) continue;
                                    if (!FindEfficiencyCut(hSigSrc, 0.90, cut90, err90)) continue;
                                    
                                    xCut.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                                    exCut.push_back(0.0);
                                    
                                    y70.push_back(cut70);
                                    ey70.push_back(err70);
                                    y80.push_back(cut80);
                                    ey80.push_back(err80);
                                    y90.push_back(cut90);
                                    ey90.push_back(err90);
                                    
                                    yMinEff = std::min(yMinEff, std::min(cut70 - err70, std::min(cut80 - err80, cut90 - err90)));
                                    yMaxEff = std::max(yMaxEff, std::max(cut70 + err70, std::max(cut80 + err80, cut90 + err90)));
                                }
                                
                                if (xCut.empty()) continue;
                                
                                if (!std::isfinite(yMinEff) || !std::isfinite(yMaxEff))
                                {
                                    yMinEff = 0.0;
                                    yMaxEff = 1.0;
                                }
                                const double padEff = (yMaxEff > yMinEff) ? (0.48 * (yMaxEff - yMinEff)) : 0.25;
                                
                                TCanvas cEff(
                                             TString::Format("c_ppg12IsoCutFits_%s_%s",
                                                             trigAA.c_str(), cb.folder.c_str()).Data(),
                                             "c_ppg12IsoCutFits", 1100, 800
                                             );
                                ApplyCanvasMargins1D(cEff);
                                cEff.SetTopMargin(0.10);
                                cEff.cd();
                                
                                TH1F hFrameEff(
                                               TString::Format("hFrame_ppg12IsoCutFits_%s_%s",
                                                               trigAA.c_str(), cb.folder.c_str()).Data(),
                                               "", 100, kPtEdges[1], kPtEdges.back()
                                               );
                                hFrameEff.SetDirectory(nullptr);
                                hFrameEff.SetStats(0);
                                hFrameEff.SetMinimum(std::max(0.0, yMinEff - padEff));
                                hFrameEff.SetMaximum(yMaxEff + padEff);
                                hFrameEff.GetXaxis()->SetTitle("Cluster p_{T} [GeV]");
                                hFrameEff.GetYaxis()->SetTitle("E_{T}^{iso} Cutoff [GeV]");
                                hFrameEff.GetXaxis()->SetTitleSize(0.060);
                                hFrameEff.GetYaxis()->SetTitleSize(0.060);
                                hFrameEff.GetXaxis()->SetLabelSize(0.050);
                                hFrameEff.GetYaxis()->SetLabelSize(0.050);
                                hFrameEff.GetYaxis()->SetTitleOffset(1.05);
                                hFrameEff.Draw();
                                
                                TGraphErrors g90((int)xCut.size(), &xCut[0], &y90[0], &exCut[0], &ey90[0]);
                                g90.SetLineWidth(2);
                                g90.SetLineColor(kMagenta + 1);
                                g90.SetMarkerColor(kMagenta + 1);
                                g90.SetMarkerStyle(20);
                                g90.SetMarkerSize(1.4);
                                
                                TGraphErrors g80((int)xCut.size(), &xCut[0], &y80[0], &exCut[0], &ey80[0]);
                                g80.SetLineWidth(2);
                                g80.SetLineColor(kGreen + 2);
                                g80.SetMarkerColor(kGreen + 2);
                                g80.SetMarkerStyle(21);
                                g80.SetMarkerSize(1.4);
                                
                                TGraphErrors g70((int)xCut.size(), &xCut[0], &y70[0], &exCut[0], &ey70[0]);
                                g70.SetLineWidth(2);
                                g70.SetLineColor(kBlue + 1);
                                g70.SetMarkerColor(kBlue + 1);
                                g70.SetMarkerStyle(22);
                                g70.SetMarkerSize(1.5);
                                
                                const double fitXLo = 14.0;
                                const double fitXHi = kPtEdges.back();
                                
                                struct IsoCutCentSummaryState
                                {
                                    vector<double> x90, ex90, y90, ey90;
                                    vector<double> x80, ex80, y80, ey80;
                                    vector<double> x70, ex70, y70, ey70;
                                };
                                static map<string, IsoCutCentSummaryState> isoCutCentSummaryStates;
                                
                                const string effSummaryKey =
                                trigAA + "_" + embeddedLabel + "_" + baseVariantPtSummaryDir;
                                auto& centSummary = isoCutCentSummaryStates[effSummaryKey];
                                if (cb.folder == centBins.front().folder)
                                {
                                    centSummary = IsoCutCentSummaryState{};
                                }
                                
                                auto ComputeFlatCutoff = [&](TGraphErrors& g, double& cutoff) -> bool
                                {
                                    cutoff = 0.0;
                                    
                                    double sumW = 0.0;
                                    double sumWY = 0.0;
                                    int nUsed = 0;
                                    
                                    for (int ip = 0; ip < g.GetN(); ++ip)
                                    {
                                        double x = 0.0, y = 0.0;
                                        g.GetPoint(ip, x, y);
                                        if (x < fitXLo || x > fitXHi) continue;
                                        
                                        const double ey = g.GetErrorY(ip);
                                        const double w = (ey > 0.0) ? (1.0 / (ey * ey)) : 1.0;
                                        sumW += w;
                                        sumWY += w * y;
                                        ++nUsed;
                                    }
                                    
                                    if (nUsed <= 0 || !(sumW > 0.0)) return false;
                                    
                                    cutoff = sumWY / sumW;
                                    return std::isfinite(cutoff);
                                };
                                
                                double flat90 = 0.0;
                                double flat80 = 0.0;
                                double flat70 = 0.0;
                                
                                const bool haveF90 = ComputeFlatCutoff(g90, flat90);
                                const bool haveF80 = ComputeFlatCutoff(g80, flat80);
                                const bool haveF70 = ComputeFlatCutoff(g70, flat70);
                                
                                const double centCenter = 0.5 * (cb.lo + cb.hi);
                                const double centHalfWidth = 0.5 * (cb.hi - cb.lo);
                                if (haveF90)
                                {
                                    centSummary.x90.push_back(centCenter);
                                    centSummary.ex90.push_back(centHalfWidth);
                                    centSummary.y90.push_back(flat90);
                                    centSummary.ey90.push_back(0.0);
                                }
                                if (haveF80)
                                {
                                    centSummary.x80.push_back(centCenter);
                                    centSummary.ex80.push_back(centHalfWidth);
                                    centSummary.y80.push_back(flat80);
                                    centSummary.ey80.push_back(0.0);
                                }
                                if (haveF70)
                                {
                                    centSummary.x70.push_back(centCenter);
                                    centSummary.ex70.push_back(centHalfWidth);
                                    centSummary.y70.push_back(flat70);
                                    centSummary.ey70.push_back(0.0);
                                }
                                
                                TF1 f90(TString::Format("f_ppg12IsoCut90_%s_%s", trigAA.c_str(), cb.folder.c_str()).Data(), "[0]", fitXLo, fitXHi);
                                TF1 f80(TString::Format("f_ppg12IsoCut80_%s_%s", trigAA.c_str(), cb.folder.c_str()).Data(), "[0]", fitXLo, fitXHi);
                                TF1 f70(TString::Format("f_ppg12IsoCut70_%s_%s", trigAA.c_str(), cb.folder.c_str()).Data(), "[0]", fitXLo, fitXHi);
                                if (haveF90) f90.SetParameter(0, flat90);
                                if (haveF80) f80.SetParameter(0, flat80);
                                if (haveF70) f70.SetParameter(0, flat70);
                                f90.SetLineColor(kMagenta + 1); f90.SetLineWidth(3); f90.SetLineStyle(2);
                                f80.SetLineColor(kGreen + 2);   f80.SetLineWidth(3); f80.SetLineStyle(2);
                                f70.SetLineColor(kBlue + 1);    f70.SetLineWidth(3); f70.SetLineStyle(2);
                                
                                // -- bump y-range so top-left labels don't collide with high-pT points --
                                {
                                    const double noPPpad = (yMaxEff > yMinEff) ? (0.75 * (yMaxEff - yMinEff)) : 0.25;
                                    hFrameEff.SetMaximum(yMaxEff + noPPpad);
                                }
                                
                                g90.Draw("PE1 SAME");
                                g80.Draw("PE1 SAME");
                                g70.Draw("PE1 SAME");
                                if (haveF90) f90.Draw("SAME");
                                if (haveF80) f80.Draw("SAME");
                                if (haveF70) f70.Draw("SAME");
                                
                                // -- legend + labels for embedded-only version (no pp overlay) --
                                {
                                    TLegend legNoPP(0.44, 0.60, 0.92, 0.79);
                                    legNoPP.SetBorderSize(0); legNoPP.SetFillStyle(0);
                                    legNoPP.SetTextFont(42); legNoPP.SetTextSize(0.030);
                                    legNoPP.AddEntry(&g90, haveF90 ? TString::Format("90%% Efficiency, E_{T}^{iso} #sim %.2f", flat90).Data() : "90% Efficiency", "ep");
                                    legNoPP.AddEntry(&g80, haveF80 ? TString::Format("80%% Efficiency, E_{T}^{iso} #sim %.2f", flat80).Data() : "80% Efficiency", "ep");
                                    legNoPP.AddEntry(&g70, haveF70 ? TString::Format("70%% Efficiency, E_{T}^{iso} #sim %.2f", flat70).Data() : "70% Efficiency", "ep");
                                    legNoPP.Draw();
                                    
                                    TLatex tInfoNoPP;
                                    tInfoNoPP.SetNDC(true); tInfoNoPP.SetTextFont(42);
                                    tInfoNoPP.SetTextAlign(13); tInfoNoPP.SetTextSize(0.030);
                                    tInfoNoPP.DrawLatex(0.18, 0.88, embeddedLabel.c_str());
                                    tInfoNoPP.DrawLatex(0.18, 0.83, TString::Format("%d-%d%% centrality", cb.lo, cb.hi).Data());
                                    tInfoNoPP.DrawLatex(0.18, 0.78,
                                                        TString::Format("|v_{z}| < %d cm,  |#eta^{#gamma}| < %.1f", kAA_VzCut, kPhotonEtaAbsMax).Data());
                                    tInfoNoPP.DrawLatex(0.18, 0.73,
                                                        TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                                    
                                    TLatex tSphNoPP;
                                    tSphNoPP.SetNDC(true); tSphNoPP.SetTextFont(42); tSphNoPP.SetTextAlign(33);
                                    tSphNoPP.SetTextSize(0.042);
                                    tSphNoPP.DrawLatex(0.92, 0.88, "#bf{sPHENIX} #it{Internal}");
                                    tSphNoPP.SetTextSize(0.034);
                                    tSphNoPP.DrawLatex(0.92, 0.83, "Pythia Emb  #sqrt{s_{NN}} = 200 GeV");
                                    
                                    // -- save embedded-only version (no pp overlay) --
                                    SaveCanvas(cEff, JoinPath(baseVariantPtSummaryDir,
                                                              TString::Format("ppg12Style_isoCutEfficiencyFits_%s_noPP.png", cb.folder.c_str()).Data()));
                                    
                                    const bool isLastCentBin = (cb.folder == centBins.back().folder);
                                    if (isLastCentBin)
                                    {
                                        double yMinCent = std::numeric_limits<double>::max();
                                        double yMaxCent = -std::numeric_limits<double>::max();
                                        auto UpdateCentRange = [&](const vector<double>& yVals, const vector<double>& eyVals)
                                        {
                                            for (std::size_t i = 0; i < yVals.size(); ++i)
                                            {
                                                yMinCent = std::min(yMinCent, yVals[i] - eyVals[i]);
                                                yMaxCent = std::max(yMaxCent, yVals[i] + eyVals[i]);
                                            }
                                        };
                                        UpdateCentRange(centSummary.y90, centSummary.ey90);
                                        UpdateCentRange(centSummary.y80, centSummary.ey80);
                                        UpdateCentRange(centSummary.y70, centSummary.ey70);
                                        
                                        if (std::isfinite(yMinCent) && std::isfinite(yMaxCent))
                                        {
                                            const double centPad = (yMaxCent > yMinCent) ? (0.65 * (yMaxCent - yMinCent)) : 0.60;
                                            
                                            TCanvas cEffCentSummary(
                                                                    TString::Format("c_ppg12IsoCutEffSummaryVsCent_%s", trigAA.c_str()).Data(),
                                                                    "c_ppg12IsoCutEffSummaryVsCent", 900, 700
                                                                    );
                                            ApplyCanvasMargins1D(cEffCentSummary);
                                            cEffCentSummary.cd();
                                            
                                            TH1F hFrameEffCent(
                                                               TString::Format("hFrame_ppg12IsoCutEffSummaryVsCent_%s", trigAA.c_str()).Data(),
                                                               "", 100, centBins.front().lo, centBins.back().hi
                                                               );
                                            hFrameEffCent.SetDirectory(nullptr);
                                            hFrameEffCent.SetStats(0);
                                            hFrameEffCent.SetMinimum(std::max(0.0, yMinCent - centPad));
                                            hFrameEffCent.SetMaximum(yMaxCent + centPad);
                                            hFrameEffCent.GetXaxis()->SetTitle("Centrality [%]");
                                            hFrameEffCent.GetYaxis()->SetTitle("E_{T}^{iso} Cutoff [GeV]");
                                            hFrameEffCent.GetXaxis()->SetTitleSize(0.055);
                                            hFrameEffCent.GetYaxis()->SetTitleSize(0.055);
                                            hFrameEffCent.GetXaxis()->SetLabelSize(0.045);
                                            hFrameEffCent.GetYaxis()->SetLabelSize(0.045);
                                            hFrameEffCent.GetYaxis()->SetTitleOffset(1.15);
                                            hFrameEffCent.Draw();
                                            
                                            TLegend legEffCent(0.48, 0.62, 0.92, 0.78);
                                            legEffCent.SetBorderSize(0);
                                            legEffCent.SetFillStyle(0);
                                            legEffCent.SetTextFont(42);
                                            legEffCent.SetTextSize(0.028);
                                            
                                            if (!centSummary.x90.empty())
                                            {
                                                std::vector<double> zeroEx90(centSummary.x90.size(), 0.0);
                                                TGraphErrors g90Cent(
                                                                     (int)centSummary.x90.size(),
                                                                     &centSummary.x90[0], &centSummary.y90[0],
                                                                     &zeroEx90[0], &centSummary.ey90[0]
                                                                     );
                                                g90Cent.SetLineWidth(2);
                                                g90Cent.SetLineColor(kMagenta + 1);
                                                g90Cent.SetMarkerColor(kMagenta + 1);
                                                g90Cent.SetMarkerStyle(20);
                                                g90Cent.SetMarkerSize(1.2);
                                                g90Cent.Draw("P SAME");
                                                TF1 *fit90 = new TF1("fit90", "pol1", centBins.front().lo, centBins.back().hi);
                                                fit90->SetLineColor(kMagenta + 1);
                                                fit90->SetLineWidth(2);
                                                fit90->SetLineStyle(2);
                                                g90Cent.Fit(fit90, "QNR");
                                                fit90->Draw("SAME");
                                                legEffCent.AddEntry(&g90Cent,
                                                                    TString::Format("90%% Eff: y = %.4fx %+.2f", fit90->GetParameter(1), fit90->GetParameter(0)).Data(), "lp");
                                                
                                                if (!centSummary.x80.empty())
                                                {
                                                    std::vector<double> zeroEx80(centSummary.x80.size(), 0.0);
                                                    TGraphErrors g80Cent(
                                                                         (int)centSummary.x80.size(),
                                                                         &centSummary.x80[0], &centSummary.y80[0],
                                                                         &zeroEx80[0], &centSummary.ey80[0]
                                                                         );
                                                    g80Cent.SetLineWidth(2);
                                                    g80Cent.SetLineColor(kGreen + 2);
                                                    g80Cent.SetMarkerColor(kGreen + 2);
                                                    g80Cent.SetMarkerStyle(21);
                                                    g80Cent.SetMarkerSize(1.2);
                                                    g80Cent.Draw("P SAME");
                                                    TF1 *fit80 = new TF1("fit80", "pol1", centBins.front().lo, centBins.back().hi);
                                                    fit80->SetLineColor(kGreen + 2);
                                                    fit80->SetLineWidth(2);
                                                    fit80->SetLineStyle(2);
                                                    g80Cent.Fit(fit80, "QNR");
                                                    fit80->Draw("SAME");
                                                    legEffCent.AddEntry(&g80Cent,
                                                                        TString::Format("80%% Eff: y = %.4fx %+.2f", fit80->GetParameter(1), fit80->GetParameter(0)).Data(), "lp");
                                                    
                                                    if (!centSummary.x70.empty())
                                                    {
                                                        std::vector<double> zeroEx70(centSummary.x70.size(), 0.0);
                                                        TGraphErrors g70Cent(
                                                                             (int)centSummary.x70.size(),
                                                                             &centSummary.x70[0], &centSummary.y70[0],
                                                                             &zeroEx70[0], &centSummary.ey70[0]
                                                                             );
                                                        g70Cent.SetLineWidth(2);
                                                        g70Cent.SetLineColor(kBlue + 1);
                                                        g70Cent.SetMarkerColor(kBlue + 1);
                                                        g70Cent.SetMarkerStyle(22);
                                                        g70Cent.SetMarkerSize(1.5);
                                                        g70Cent.Draw("P SAME");
                                                        TF1 *fit70 = new TF1("fit70", "pol1", centBins.front().lo, centBins.back().hi);
                                                        fit70->SetLineColor(kBlue + 1);
                                                        fit70->SetLineWidth(2);
                                                        fit70->SetLineStyle(2);
                                                        g70Cent.Fit(fit70, "QNR");
                                                        fit70->Draw("SAME");
                                                        legEffCent.AddEntry(&g70Cent,
                                                                            TString::Format("70%% Eff: y = %.4fx %+.2f", fit70->GetParameter(1), fit70->GetParameter(0)).Data(), "lp");
                                                        
                                                        legEffCent.Draw();
                                                        
                                                        TLatex tTitleCent;
                                                        tTitleCent.SetNDC(true);
                                                        tTitleCent.SetTextFont(42);
                                                        tTitleCent.SetTextAlign(23);
                                                        tTitleCent.SetTextSize(0.042);
                                                        tTitleCent.DrawLatex(0.50, 0.98, "Flat E_{T}^{iso} cutoff vs centrality");
                                                        
                                                        TLatex tInfoCent;
                                                        tInfoCent.SetNDC(true);
                                                        tInfoCent.SetTextFont(42);
                                                        tInfoCent.SetTextAlign(13);
                                                        tInfoCent.SetTextSize(0.030);
                                                        tInfoCent.DrawLatex(0.18, 0.88, embeddedLabel.c_str());
                                                        tInfoCent.DrawLatex(0.18, 0.83,
                                                                            TString::Format("linear fit over p_{T}^{#gamma} #geq %.0f GeV", fitXLo).Data());
                                                        tInfoCent.DrawLatex(0.18, 0.78,
                                                                            TString::Format("|v_{z}| < %d cm,  |#eta^{#gamma}| < %.1f", kAA_VzCut, kPhotonEtaAbsMax).Data());
                                                        tInfoCent.DrawLatex(0.18, 0.73,
                                                                            TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                                                        
                                                        TLatex tSphCent;
                                                        tSphCent.SetNDC(true);
                                                        tSphCent.SetTextFont(42);
                                                        tSphCent.SetTextAlign(33);
                                                        tSphCent.SetTextSize(0.042);
                                                        tSphCent.DrawLatex(0.92, 0.88, "#bf{sPHENIX} #it{Internal}");
                                                        tSphCent.SetTextSize(0.034);
                                                        tSphCent.DrawLatex(0.92, 0.83, "Pythia Emb  #sqrt{s_{NN}} = 200 GeV");
                                                        
                                                        SaveCanvas(cEffCentSummary, JoinPath(baseVariantPtSummaryDir,
                                                                                             "ppg12Style_isoCutEfficiencyFits_vsCentrality_noPP.png"));
                                                    }
                                                    else
                                                    {
                                                        legEffCent.Draw();
                                                        TLatex tTitleCent;
                                                        tTitleCent.SetNDC(true);
                                                        tTitleCent.SetTextFont(42);
                                                        tTitleCent.SetTextAlign(23);
                                                        tTitleCent.SetTextSize(0.042);
                                                        tTitleCent.DrawLatex(0.50, 0.98, "Flat E_{T}^{iso} cutoff vs centrality");
                                                        SaveCanvas(cEffCentSummary, JoinPath(baseVariantPtSummaryDir,
                                                                                             "ppg12Style_isoCutEfficiencyFits_vsCentrality_noPP.png"));
                                                    }
                                                }
                                                else
                                                {
                                                    legEffCent.Draw();
                                                    TLatex tTitleCent;
                                                    tTitleCent.SetNDC(true);
                                                    tTitleCent.SetTextFont(42);
                                                    tTitleCent.SetTextAlign(23);
                                                    tTitleCent.SetTextSize(0.042);
                                                    tTitleCent.DrawLatex(0.50, 0.98, "Flat E_{T}^{iso} cutoff vs centrality");
                                                    SaveCanvas(cEffCentSummary, JoinPath(baseVariantPtSummaryDir,
                                                                                         "ppg12Style_isoCutEfficiencyFits_vsCentrality_noPP.png"));
                                                }
                                            }
                                            
                                            isoCutCentSummaryStates.erase(effSummaryKey);
                                        }
                                    }
                                }
                                
                                // -- pp reference (matched pp SIM, open markers) --
                                vector<double> xCutPP, exCutPP, y90PP, ey90PP, y80PP, ey80PP, y70PP, ey70PP;
                                if (matchedPPSimTop)
                                {
                                    for (int iptPP = 1; iptPP < kNPtBins; ++iptPP)
                                    {
                                        const PtBin& bPP = PtBins()[iptPP];
                                        const string hPPSigName = "h_EisoReco_truthSigMatched" + bPP.suffix;
                                        TH1* hPPSig = dynamic_cast<TH1*>(matchedPPSimTop->Get(hPPSigName.c_str()));
                                        if (!hPPSig || hPPSig->GetEntries() <= 0.0) continue;
                                        double c70=0,e70=0, c80=0,e80=0, c90=0,e90=0;
                                        if (!FindEfficiencyCut(hPPSig, 0.70, c70, e70)) continue;
                                        if (!FindEfficiencyCut(hPPSig, 0.80, c80, e80)) continue;
                                        if (!FindEfficiencyCut(hPPSig, 0.90, c90, e90)) continue;
                                        xCutPP.push_back(0.5*(kPtEdges[(std::size_t)iptPP]+kPtEdges[(std::size_t)iptPP+1]));
                                        exCutPP.push_back(0.0);
                                        y90PP.push_back(c90); ey90PP.push_back(e90);
                                        y80PP.push_back(c80); ey80PP.push_back(e80);
                                        y70PP.push_back(c70); ey70PP.push_back(e70);
                                    }
                                }
                                // update y-range to include pp points
                                for (std::size_t ipp = 0; ipp < xCutPP.size(); ++ipp)
                                {
                                    yMinEff = std::min(yMinEff, std::min(y70PP[ipp]-ey70PP[ipp], std::min(y80PP[ipp]-ey80PP[ipp], y90PP[ipp]-ey90PP[ipp])));
                                    yMaxEff = std::max(yMaxEff, std::max(y70PP[ipp]+ey70PP[ipp], std::max(y80PP[ipp]+ey80PP[ipp], y90PP[ipp]+ey90PP[ipp])));
                                }
                                if (!xCutPP.empty())
                                {
                                    const double padEffUpd = (yMaxEff > yMinEff) ? (0.48 * (yMaxEff - yMinEff)) : 0.25;
                                    hFrameEff.SetMinimum(std::max(0.0, yMinEff - padEffUpd));
                                    hFrameEff.SetMaximum(yMaxEff + padEffUpd);
                                }
                                TGraphErrors *gPP90=nullptr, *gPP80=nullptr, *gPP70=nullptr;
                                if (!xCutPP.empty())
                                {
                                    gPP90 = new TGraphErrors((int)xCutPP.size(), &xCutPP[0], &y90PP[0], &exCutPP[0], &ey90PP[0]);
                                    gPP90->SetLineWidth(2); gPP90->SetLineColor(kMagenta+1); gPP90->SetMarkerColor(kMagenta+1);
                                    gPP90->SetMarkerStyle(24); gPP90->SetMarkerSize(1.4);
                                    gPP90->Draw("PE1 SAME");
                                    gPP80 = new TGraphErrors((int)xCutPP.size(), &xCutPP[0], &y80PP[0], &exCutPP[0], &ey80PP[0]);
                                    gPP80->SetLineWidth(2); gPP80->SetLineColor(kGreen+2); gPP80->SetMarkerColor(kGreen+2);
                                    gPP80->SetMarkerStyle(25); gPP80->SetMarkerSize(1.4);
                                    gPP80->Draw("PE1 SAME");
                                    gPP70 = new TGraphErrors((int)xCutPP.size(), &xCutPP[0], &y70PP[0], &exCutPP[0], &ey70PP[0]);
                                    gPP70->SetLineWidth(2); gPP70->SetLineColor(kBlue+1); gPP70->SetMarkerColor(kBlue+1);
                                    gPP70->SetMarkerStyle(26); gPP70->SetMarkerSize(1.5);
                                    gPP70->Draw("PE1 SAME");
                                }
                                
                                TLegend legEff(0.38, 0.58, 0.92, 0.78);
                                legEff.SetBorderSize(0);
                                legEff.SetFillStyle(0);
                                legEff.SetTextFont(42);
                                legEff.SetTextSize(0.024);
                                legEff.SetNColumns(2);
                                legEff.SetColumnSeparation(-0.02);
                                legEff.AddEntry(&g90, TString::Format("90%% Efficiency (%s)", embLegTag.c_str()).Data(), "ep");
                                if (gPP90) legEff.AddEntry(gPP90, TString::Format("90%% Efficiency (%s)", matchedPPLabel.c_str()).Data(), "ep");
                                legEff.AddEntry(&g80, TString::Format("80%% Efficiency (%s)", embLegTag.c_str()).Data(), "ep");
                                if (gPP80) legEff.AddEntry(gPP80, TString::Format("80%% Efficiency (%s)", matchedPPLabel.c_str()).Data(), "ep");
                                legEff.AddEntry(&g70, TString::Format("70%% Efficiency (%s)", embLegTag.c_str()).Data(), "ep");
                                if (gPP70) legEff.AddEntry(gPP70, TString::Format("70%% Efficiency (%s)", matchedPPLabel.c_str()).Data(), "ep");
                                legEff.Draw();
                                
                                TLatex tInfoEff;
                                tInfoEff.SetNDC(true);
                                tInfoEff.SetTextFont(42);
                                tInfoEff.SetTextAlign(13);
                                tInfoEff.SetTextSize(0.030);
                                tInfoEff.DrawLatex(0.18, 0.88, embeddedLabel.c_str());
                                tInfoEff.DrawLatex(0.18, 0.83, TString::Format("%d-%d%% centrality", cb.lo, cb.hi).Data());
                                tInfoEff.DrawLatex(0.18, 0.78,
                                                   TString::Format("|v_{z}| < %d cm,  |#eta^{#gamma}| < %.1f", kAA_VzCut, kPhotonEtaAbsMax).Data());
                                
                                tInfoEff.DrawLatex(0.18, 0.73,
                                                   TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                                
                                {
                                    TLatex tSph;
                                    tSph.SetNDC(true);
                                    tSph.SetTextFont(42);
                                    tSph.SetTextAlign(33);
                                    tSph.SetTextSize(0.042);
                                    tSph.DrawLatex(0.92, 0.88, "#bf{sPHENIX} #it{Internal}");
                                    tSph.SetTextSize(0.034);
                                    tSph.DrawLatex(0.92, 0.83, "Pythia Emb  #sqrt{s_{NN}} = 200 GeV");
                                }
                                
                                SaveCanvas(cEff, JoinPath(baseVariantPtSummaryDir,
                                                          TString::Format("ppg12Style_isoCutEfficiencyFits_%s.png", cb.folder.c_str()).Data()));
                                if (gPP90) delete gPP90;
                                if (gPP80) delete gPP80;
                                if (gPP70) delete gPP70;
                            }
                            
                            if (fMatchedPPSim) { fMatchedPPSim->Close(); delete fMatchedPPSim; fMatchedPPSim = nullptr; }
                            
                            // -- standalone pp 5+10+20 merged pythia-only iso-cut efficiency plot (no centrality) --
                            if (ppSimMergedTop)
                            {
                                vector<double> xPPM, exPPM, y90M, ey90M, y80M, ey80M, y70M, ey70M;
                                double yMinM = std::numeric_limits<double>::max();
                                double yMaxM = -std::numeric_limits<double>::max();
                                for (int ipt = 1; ipt < kNPtBins; ++ipt)
                                {
                                    const PtBin& b = PtBins()[ipt];
                                    const string hName = "h_EisoReco_truthSigMatched" + b.suffix;
                                    TH1* hSig = dynamic_cast<TH1*>(ppSimMergedTop->Get(hName.c_str()));
                                    if (!hSig || hSig->GetEntries() <= 0.0) continue;
                                    double c70=0,e70=0, c80=0,e80=0, c90=0,e90=0;
                                    if (!FindEfficiencyCut(hSig, 0.70, c70, e70)) continue;
                                    if (!FindEfficiencyCut(hSig, 0.80, c80, e80)) continue;
                                    if (!FindEfficiencyCut(hSig, 0.90, c90, e90)) continue;
                                    xPPM.push_back(0.5*(kPtEdges[(std::size_t)ipt]+kPtEdges[(std::size_t)ipt+1]));
                                    exPPM.push_back(0.0);
                                    y90M.push_back(c90); ey90M.push_back(e90);
                                    y80M.push_back(c80); ey80M.push_back(e80);
                                    y70M.push_back(c70); ey70M.push_back(e70);
                                    yMinM = std::min(yMinM, std::min(c70-e70, std::min(c80-e80, c90-e90)));
                                    yMaxM = std::max(yMaxM, std::max(c70+e70, std::max(c80+e80, c90+e90)));
                                }
                                if (!xPPM.empty())
                                {
                                    if (!std::isfinite(yMinM) || !std::isfinite(yMaxM)) { yMinM = 0.0; yMaxM = 1.0; }
                                    const double padM = (yMaxM > yMinM) ? (1.2 * (yMaxM - yMinM)) : 0.25;
                                    TCanvas cPPM("c_ppg12IsoCutFits_ppMerged", "", 900, 700);
                                    ApplyCanvasMargins1D(cPPM);
                                    cPPM.SetTopMargin(0.10);
                                    cPPM.cd();
                                    TH1F hFrM("hFrM_ppMerged", "", 100, kPtEdges[1], kPtEdges[kNPtBins]);
                                    hFrM.SetDirectory(nullptr); hFrM.SetStats(0);
                                    hFrM.SetMinimum(std::max(0.0, yMinM - padM));
                                    hFrM.SetMaximum(yMaxM + padM);
                                    hFrM.GetXaxis()->SetTitle("Cluster p_{T} [GeV]");
                                    hFrM.GetYaxis()->SetTitle("E_{T}^{iso} Cutoff [GeV]");
                                    hFrM.GetXaxis()->SetTitleSize(0.060);
                                    hFrM.GetYaxis()->SetTitleSize(0.060);
                                    hFrM.GetXaxis()->SetLabelSize(0.050);
                                    hFrM.GetYaxis()->SetLabelSize(0.050);
                                    hFrM.GetYaxis()->SetTitleOffset(1.05);
                                    hFrM.Draw();
                                    
                                    TGraphErrors gM90((int)xPPM.size(), &xPPM[0], &y90M[0], &exPPM[0], &ey90M[0]);
                                    gM90.SetLineWidth(2); gM90.SetLineColor(kMagenta+1); gM90.SetMarkerColor(kMagenta+1);
                                    gM90.SetMarkerStyle(20); gM90.SetMarkerSize(1.4); gM90.Draw("PE1 SAME");
                                    TGraphErrors gM80((int)xPPM.size(), &xPPM[0], &y80M[0], &exPPM[0], &ey80M[0]);
                                    gM80.SetLineWidth(2); gM80.SetLineColor(kGreen+2); gM80.SetMarkerColor(kGreen+2);
                                    gM80.SetMarkerStyle(21); gM80.SetMarkerSize(1.4); gM80.Draw("PE1 SAME");
                                    TGraphErrors gM70((int)xPPM.size(), &xPPM[0], &y70M[0], &exPPM[0], &ey70M[0]);
                                    gM70.SetLineWidth(2); gM70.SetLineColor(kBlue+1); gM70.SetMarkerColor(kBlue+1);
                                    gM70.SetMarkerStyle(22); gM70.SetMarkerSize(1.5); gM70.Draw("PE1 SAME");
                                    
                                    const double fitXLoM = kPtEdges[1];
                                    const double fitXHiM = kPtEdges.back();
                                    TF1 fM90("fM90_ppMerged", "pol1", fitXLoM, fitXHiM);
                                    TF1 fM80("fM80_ppMerged", "pol1", fitXLoM, fitXHiM);
                                    TF1 fM70("fM70_ppMerged", "pol1", fitXLoM, fitXHiM);
                                    fM90.SetLineColor(kMagenta+1); fM90.SetLineWidth(3);
                                    fM80.SetLineColor(kGreen+2);   fM80.SetLineWidth(3);
                                    fM70.SetLineColor(kBlue+1);    fM70.SetLineWidth(3);
                                    if (gM90.GetN() >= 2) { gM90.Fit(&fM90, "Q0"); fM90.Draw("SAME"); }
                                    if (gM80.GetN() >= 2) { gM80.Fit(&fM80, "Q0"); fM80.Draw("SAME"); }
                                    if (gM70.GetN() >= 2) { gM70.Fit(&fM70, "Q0"); fM70.Draw("SAME"); }
                                    
                                    TLegend legM(0.50, 0.62, 0.92, 0.78);
                                    legM.SetBorderSize(0); legM.SetFillStyle(0);
                                    legM.SetTextFont(42); legM.SetTextSize(0.030);
                                    legM.AddEntry(&gM90, "90% Efficiency", "ep");
                                    legM.AddEntry(&gM80, "80% Efficiency", "ep");
                                    legM.AddEntry(&gM70, "70% Efficiency", "ep");
                                    legM.Draw();
                                    
                                    TLatex tInfoM;
                                    tInfoM.SetNDC(true); tInfoM.SetTextFont(42);
                                    tInfoM.SetTextAlign(13); tInfoM.SetTextSize(0.030);
                                    tInfoM.DrawLatex(0.18, 0.88, "Pythia,  #sqrt{s} = 200 GeV");
                                    if (gM90.GetN() >= 2)
                                        tInfoM.DrawLatex(0.18, 0.83,
                                                         TString::Format("90%%: E_{T}^{iso} = %.3f  + %.3fp_{T}", fM90.GetParameter(0), fM90.GetParameter(1)).Data());
                                    if (gM80.GetN() >= 2)
                                        tInfoM.DrawLatex(0.18, 0.78,
                                                         TString::Format("80%%: E_{T}^{iso} = %.3f  + %.3fp_{T}", fM80.GetParameter(0), fM80.GetParameter(1)).Data());
                                    if (gM70.GetN() >= 2)
                                        tInfoM.DrawLatex(0.18, 0.73,
                                                         TString::Format("70%%: E_{T}^{iso} = %.3f  + %.3fp_{T}", fM70.GetParameter(0), fM70.GetParameter(1)).Data());
                                    tInfoM.DrawLatex(0.18, 0.68,
                                                     TString::Format("vtx |z| < %d cm, |#eta^{#gamma}|<%.1f", kAA_VzCut, kPhotonEtaAbsMax).Data());
                                    
                                    {
                                        TLatex tSph;
                                        tSph.SetNDC(true); tSph.SetTextFont(42); tSph.SetTextAlign(33);
                                        tSph.SetTextSize(0.042);
                                        tSph.DrawLatex(0.92, 0.88, "#bf{sPHENIX} #it{Internal}");
                                    }
                                    
                                    SaveCanvas(cPPM, JoinPath(baseVariantPtSummaryDir,
                                                              "ppg12Style_isoCutEfficiencyFits_ppMerged.png"));
                                }
                            }
                            
                            fPhoCut->Close();
                            
                            delete fPhoCut;
                            fPhoCut = nullptr;
                        }
                        else
                        {
                            if (fPhoCut) { fPhoCut->Close(); delete fPhoCut; fPhoCut = nullptr; }
                            cout << ANSI_BOLD_YEL
                            << "[WARN] Missing photon+jet embedded baseVariant input for PPG12-style iso-cut fits: "
                            << phoCutIn
                            << ANSI_RESET << "\n";
                        }
                    }
                    else
                    {
                        cout << ANSI_BOLD_YEL
                        << "[WARN] Photon+jet embedded SIM toggle/path is unavailable for baseVariant PPG12-style iso-cut fits."
                        << ANSI_RESET << "\n";
                    }
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
                            {
                                const string variantGaussPtDir = JoinPath(JoinPath(ueCompBase, H.variant), "ppAuAu_meanIsoEt_versusPT");
                                EnsureDir(variantGaussPtDir);
                                SaveCanvas(cGM, JoinPath(variantGaussPtDir,
                                                         TString::Format("gaussMean_auauCentOverlay_vs_pT_%s.png", H.variant.c_str()).Data()));
                            }
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
        } // end if (!perVariantIsoQAOnlyActive) — pT/cent summary sections
        
        if (!forEmbeddedSim)
        {
#if 0
            {
                const int centPtColors[] = {kBlack, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 1,
                                            kCyan + 1, kYellow + 2, kViolet + 1};
                
                // --- Gaussian mean vs pT: all centralities overlaid + pp ---
                {
                    vector<TGraphErrors*> gCentMean;
                    vector<string>        centMeanLabels;
                    vector<double>        xPPMean, exPPMean, yPPMean, eyPPMean;
                    double yMinMean = std::numeric_limits<double>::max();
                    double yMaxMean = -std::numeric_limits<double>::max();
                    
                    for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                    {
                        const auto& cb = centBins[ic];
                        vector<double> xP, exP, yP, eyP;
                        
                        for (int ipt = 0; ipt < kNPtBins; ++ipt)
                        {
                            if (!gaussFilled[ivH][ipt][ic]) continue;
                            
                            xP.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                            exP.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                            yP.push_back(gaussMean[ivH][ipt][ic]);
                            eyP.push_back(gaussMeanErr[ivH][ipt][ic]);
                            
                            yMinMean = std::min(yMinMean, yP.back() - eyP.back());
                            yMaxMean = std::max(yMaxMean, yP.back() + eyP.back());
                        }
                        
                        if (xP.empty()) continue;
                        
                        TGraphErrors* g = new TGraphErrors((int)xP.size(), &xP[0], &yP[0], &exP[0], &eyP[0]);
                        const int ci = (ic < 8) ? centPtColors[ic] : kBlack;
                        g->SetLineWidth(2);
                        g->SetLineColor(ci);
                        g->SetMarkerColor(ci);
                        g->SetMarkerStyle(20);
                        g->SetMarkerSize(1.2);
                        gCentMean.push_back(g);
                        centMeanLabels.push_back(TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                    }
                    
                    for (int ipt = 0; ipt < kNPtBins; ++ipt)
                    {
                        if (!ppGaussFilledByPt[ipt]) continue;
                        
                        xPPMean.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                        exPPMean.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                        yPPMean.push_back(ppGaussMeanByPt[ipt]);
                        eyPPMean.push_back(ppGaussMeanErrByPt[ipt]);
                        
                        yMinMean = std::min(yMinMean, yPPMean.back() - eyPPMean.back());
                        yMaxMean = std::max(yMaxMean, yPPMean.back() + eyPPMean.back());
                    }
                    
                    if (!gCentMean.empty() || !xPPMean.empty())
                    {
                        if (!std::isfinite(yMinMean) || !std::isfinite(yMaxMean))
                        {
                            yMinMean = 0.0;
                            yMaxMean = 1.0;
                        }
                        const double padMean = (yMaxMean > yMinMean) ? (0.45 * (yMaxMean - yMinMean)) : 0.5;
                        
                        TCanvas cMeanPtCent(
                            TString::Format("c_gaussMean_ppCentOverlay_vsPt_%s_%s",
                                            trigAA.c_str(), H.variant.c_str()).Data(),
                            "c_gaussMean_ppCentOverlay_vsPt", 900, 700);
                        ApplyCanvasMargins1D(cMeanPtCent);
                        cMeanPtCent.cd();
                        
                        TH1F hFrameMean(
                            TString::Format("hFrame_gaussMean_ppCentOverlay_vsPt_%s_%s",
                                            trigAA.c_str(), H.variant.c_str()).Data(),
                            "", 100, kPtEdges.front(), kPtEdges.back());
                        hFrameMean.SetDirectory(nullptr);
                        hFrameMean.SetStats(0);
                        hFrameMean.SetMinimum(0.0);
                        hFrameMean.SetMaximum(yMaxMean + padMean);
                        hFrameMean.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hFrameMean.GetYaxis()->SetTitle("Gaussian #mu [GeV]");
                        hFrameMean.GetXaxis()->SetTitleSize(0.055);
                        hFrameMean.GetYaxis()->SetTitleSize(0.055);
                        hFrameMean.GetXaxis()->SetLabelSize(0.045);
                        hFrameMean.GetYaxis()->SetLabelSize(0.045);
                        hFrameMean.GetYaxis()->SetTitleOffset(1.15);
                        hFrameMean.Draw();
                        
                        TGraphErrors* gPPMean = nullptr;
                        if (!xPPMean.empty())
                        {
                            gPPMean = new TGraphErrors((int)xPPMean.size(), &xPPMean[0], &yPPMean[0], &exPPMean[0], &eyPPMean[0]);
                            gPPMean->SetLineWidth(2);
                            gPPMean->SetLineColor(kRed + 1);
                            gPPMean->SetMarkerColor(kRed + 1);
                            gPPMean->SetMarkerStyle(24);
                            gPPMean->SetMarkerSize(1.3);
                            gPPMean->Draw("PE1 SAME");
                        }
                        
                        for (auto* g : gCentMean) g->Draw("PE1 SAME");
                        
                        TLegend legMean(0.18, 0.66, 0.62, 0.88);
                        legMean.SetBorderSize(0);
                        legMean.SetFillStyle(0);
                        legMean.SetTextFont(42);
                        legMean.SetTextSize(0.030);
                        legMean.SetNColumns(2);
                        if (gPPMean) legMean.AddEntry(gPPMean, "pp", "ep");
                        for (std::size_t ig = 0; ig < gCentMean.size(); ++ig)
                            legMean.AddEntry(gCentMean[ig], centMeanLabels[ig].c_str(), "ep");
                        legMean.Draw();
                        
                        TLatex tMeanTitle;
                        tMeanTitle.SetNDC(true);
                        tMeanTitle.SetTextFont(42);
                        tMeanTitle.SetTextAlign(23);
                        tMeanTitle.SetTextSize(0.042);
                        tMeanTitle.DrawLatex(0.50, 0.98,
                                             TString::Format("Gaussian #mu vs p_{T}^{#gamma}, centrality overlay + pp (%s)", H.label.c_str()).Data());
                        
                        TLatex tMeanInfo;
                        tMeanInfo.SetNDC(true);
                        tMeanInfo.SetTextFont(42);
                        tMeanInfo.SetTextAlign(13);
                        tMeanInfo.SetTextSize(0.028);
                        tMeanInfo.DrawLatex(0.18, 0.58, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
                        tMeanInfo.DrawLatex(0.18, 0.54, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                        tMeanInfo.DrawLatex(0.18, 0.50, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        tMeanInfo.DrawLatex(0.18, 0.46, TString::Format("UE: %s", H.label.c_str()).Data());
                        
                        SaveCanvas(cMeanPtCent, JoinPath(variantDir,
                                                         TString::Format("gaussMean_pp_vs_auauCentOverlay_vs_pT_%s.png", H.variant.c_str()).Data()));
                        
                        if (gPPMean) delete gPPMean;
                    }
                    
                    for (auto* g : gCentMean) delete g;
                }
                
                // --- Gaussian sigma vs pT: all centralities overlaid + pp ---
                {
                    vector<TGraphErrors*> gCentSigma;
                    vector<string>        centSigmaLabels;
                    vector<double>        xPPSigma, exPPSigma, yPPSigma, eyPPSigma;
                    double yMinSigma = std::numeric_limits<double>::max();
                    double yMaxSigma = -std::numeric_limits<double>::max();
                    
                    for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                    {
                        const auto& cb = centBins[ic];
                        vector<double> xP, exP, yP, eyP;
                        
                        for (int ipt = 0; ipt < kNPtBins; ++ipt)
                        {
                            if (!gaussFilled[ivH][ipt][ic]) continue;
                            
                            xP.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                            exP.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                            yP.push_back(gaussSigma[ivH][ipt][ic]);
                            eyP.push_back(gaussSigmaErr[ivH][ipt][ic]);
                            
                            yMinSigma = std::min(yMinSigma, yP.back() - eyP.back());
                            yMaxSigma = std::max(yMaxSigma, yP.back() + eyP.back());
                        }
                        
                        if (xP.empty()) continue;
                        
                        TGraphErrors* g = new TGraphErrors((int)xP.size(), &xP[0], &yP[0], &exP[0], &eyP[0]);
                        const int ci = (ic < 8) ? centPtColors[ic] : kBlack;
                        g->SetLineWidth(2);
                        g->SetLineColor(ci);
                        g->SetMarkerColor(ci);
                        g->SetMarkerStyle(20);
                        g->SetMarkerSize(1.2);
                        gCentSigma.push_back(g);
                        centSigmaLabels.push_back(TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                    }
                    
                    for (int ipt = 0; ipt < kNPtBins; ++ipt)
                    {
                        if (!ppGaussFilledByPt[ipt]) continue;
                        
                        xPPSigma.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                        exPPSigma.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                        yPPSigma.push_back(ppGaussSigmaByPt[ipt]);
                        eyPPSigma.push_back(ppGaussSigmaErrByPt[ipt]);
                        
                        yMinSigma = std::min(yMinSigma, yPPSigma.back() - eyPPSigma.back());
                        yMaxSigma = std::max(yMaxSigma, yPPSigma.back() + eyPPSigma.back());
                    }
                    
                    if (!gCentSigma.empty() || !xPPSigma.empty())
                    {
                        if (!std::isfinite(yMinSigma) || !std::isfinite(yMaxSigma))
                        {
                            yMinSigma = 0.0;
                            yMaxSigma = 1.0;
                        }
                        const double padSigma = (yMaxSigma > yMinSigma) ? (0.95 * (yMaxSigma - yMinSigma)) : 0.5;
                        
                        TCanvas cSigmaPtCent(
                            TString::Format("c_gaussSigma_ppCentOverlay_vsPt_%s_%s",
                                            trigAA.c_str(), H.variant.c_str()).Data(),
                            "c_gaussSigma_ppCentOverlay_vsPt", 900, 700);
                        ApplyCanvasMargins1D(cSigmaPtCent);
                        cSigmaPtCent.cd();
                        
                        TH1F hFrameSigma(
                            TString::Format("hFrame_gaussSigma_ppCentOverlay_vsPt_%s_%s",
                                            trigAA.c_str(), H.variant.c_str()).Data(),
                            "", 100, kPtEdges.front(), kPtEdges.back());
                        hFrameSigma.SetDirectory(nullptr);
                        hFrameSigma.SetStats(0);
                        hFrameSigma.SetMinimum(0.0);
                        hFrameSigma.SetMaximum(yMaxSigma + padSigma);
                        hFrameSigma.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hFrameSigma.GetYaxis()->SetTitle("Gaussian #sigma [GeV]");
                        hFrameSigma.GetXaxis()->SetTitleSize(0.055);
                        hFrameSigma.GetYaxis()->SetTitleSize(0.055);
                        hFrameSigma.GetXaxis()->SetLabelSize(0.045);
                        hFrameSigma.GetYaxis()->SetLabelSize(0.045);
                        hFrameSigma.GetYaxis()->SetTitleOffset(1.15);
                        hFrameSigma.Draw();
                        
                        TGraphErrors* gPPSigma = nullptr;
                        if (!xPPSigma.empty())
                        {
                            gPPSigma = new TGraphErrors((int)xPPSigma.size(), &xPPSigma[0], &yPPSigma[0], &exPPSigma[0], &eyPPSigma[0]);
                            gPPSigma->SetLineWidth(2);
                            gPPSigma->SetLineColor(kRed + 1);
                            gPPSigma->SetMarkerColor(kRed + 1);
                            gPPSigma->SetMarkerStyle(24);
                            gPPSigma->SetMarkerSize(1.3);
                            gPPSigma->Draw("PE1 SAME");
                        }
                        
                        for (auto* g : gCentSigma) g->Draw("PE1 SAME");
                        
                        TLegend legSigma(0.18, 0.66, 0.62, 0.88);
                        legSigma.SetBorderSize(0);
                        legSigma.SetFillStyle(0);
                        legSigma.SetTextFont(42);
                        legSigma.SetTextSize(0.030);
                        legSigma.SetNColumns(2);
                        if (gPPSigma) legSigma.AddEntry(gPPSigma, "pp", "ep");
                        for (std::size_t ig = 0; ig < gCentSigma.size(); ++ig)
                            legSigma.AddEntry(gCentSigma[ig], centSigmaLabels[ig].c_str(), "ep");
                        legSigma.Draw();
                        
                        TLatex tSigmaTitle;
                        tSigmaTitle.SetNDC(true);
                        tSigmaTitle.SetTextFont(42);
                        tSigmaTitle.SetTextAlign(23);
                        tSigmaTitle.SetTextSize(0.042);
                        tSigmaTitle.DrawLatex(0.50, 0.98,
                                              TString::Format("Gaussian #sigma vs p_{T}^{#gamma}, centrality overlay + pp (%s)", H.label.c_str()).Data());
                        
                        TLatex tSigmaInfo;
                        tSigmaInfo.SetNDC(true);
                        tSigmaInfo.SetTextFont(42);
                        tSigmaInfo.SetTextAlign(13);
                        tSigmaInfo.SetTextSize(0.028);
                        tSigmaInfo.DrawLatex(0.18, 0.58, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
                        tSigmaInfo.DrawLatex(0.18, 0.54, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                        tSigmaInfo.DrawLatex(0.18, 0.50, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        tSigmaInfo.DrawLatex(0.18, 0.46, TString::Format("UE: %s", H.label.c_str()).Data());
                        
                        SaveCanvas(cSigmaPtCent, JoinPath(variantDir,
                                                          TString::Format("gaussSigma_pp_vs_auauCentOverlay_vs_pT_%s.png", H.variant.c_str()).Data()));
                        
                        if (gPPSigma) delete gPPSigma;
                    }
                    
                    for (auto* g : gCentSigma) delete g;
                }
            }
            
#endif
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
                    const string summaryOutputDir = JoinPath(variantBaseDir, "summaryOutput");
                    const string meanSummaryDir = JoinPath(summaryOutputDir, "meanSummary");
                    const string sigmaSummaryDir = JoinPath(summaryOutputDir, "sigmaSummary");
                    EnsureDir(summaryOutputDir);
                    EnsureDir(meanSummaryDir);
                    EnsureDir(sigmaSummaryDir);
                    
                    for (const auto& cb : centBins)
                    {
                        const string centDir = JoinPath(variantBaseDir, cb.folder);
                        EnsureDir(centDir);
                        
                        // Pre-pass: accumulate Gaussian means for all 4 series across pT bins
                        vector<double> tntSubTDX, tntSubTDY, tntSubTDEY, tntSubTDSigY, tntSubTDSigEY;
                        vector<double> tntSubNTDX, tntSubNTDY, tntSubNTDEY, tntSubNTDSigY, tntSubNTDSigEY;
                        vector<double> tntSubTMX, tntSubTMY, tntSubTMEY, tntSubTMSigY, tntSubTMSigEY;
                        vector<double> tntSubNTMX, tntSubNTMY, tntSubNTMEY, tntSubNTMSigY, tntSubNTMSigEY;
                        for (int iptPre = 0; iptPre < kNPtBins; ++iptPre)
                        {
                            const PtBin& bPre = PtBins()[iptPre];
                            if (bPre.lo < 10) continue;
                            const double ptC = 0.5 * (kPtEdges[(std::size_t)iptPre] + kPtEdges[(std::size_t)iptPre + 1]);
                            const string hTN  = "h_Eiso_tight" + bPre.suffix + cb.suffix;
                            const string hNTN = "h_Eiso_nonTight" + bPre.suffix + cb.suffix;
                            auto FitAndPush = [&](TDirectory* dir, const string& hName, vector<double>& vx, vector<double>& vy, vector<double>& vey,
                                                  vector<double>& vsy, vector<double>& vsey) {
                                TH1* hSrc = dynamic_cast<TH1*>(dir->Get(hName.c_str()));
                                if (!hSrc) return;
                                TH1* hTmp = CloneTH1(hSrc, TString::Format("hPreFit_%s_%d_%p", hName.c_str(), iptPre, (void*)dir).Data());
                                if (!hTmp) return;
                                hTmp->Rebin(10); EnsureSumw2(hTmp);
                                double gM, gS, gME, gSE;
                                if (FitGaussianIterative(hTmp, gM, gS, gME, gSE))
                                { vx.push_back(ptC); vy.push_back(gM); vey.push_back(gME); vsy.push_back(gS); vsey.push_back(gSE); }
                                delete hTmp;
                            };
                            FitAndPush(dataTop, hTN,  tntSubTDX,  tntSubTDY,  tntSubTDEY,  tntSubTDSigY,  tntSubTDSigEY);
                            FitAndPush(dataTop, hNTN, tntSubNTDX, tntSubNTDY, tntSubNTDEY, tntSubNTDSigY, tntSubNTDSigEY);
                            FitAndPush(mcTop,   hTN,  tntSubTMX,  tntSubTMY,  tntSubTMEY,  tntSubTMSigY,  tntSubTMSigEY);
                            FitAndPush(mcTop,   hNTN, tntSubNTMX, tntSubNTMY, tntSubNTMEY, tntSubNTMSigY, tntSubNTMSigEY);
                        }
                        double tntSubYLo = 1e30, tntSubYHi = -1e30;
                        auto TntUpdateRange = [&](const vector<double>& y, const vector<double>& ey) {
                            for (std::size_t i = 0; i < y.size(); ++i) {
                                tntSubYLo = std::min(tntSubYLo, y[i] - ey[i]);
                                tntSubYHi = std::max(tntSubYHi, y[i] + ey[i]); } };
                        TntUpdateRange(tntSubTDY, tntSubTDEY);
                        TntUpdateRange(tntSubNTDY, tntSubNTDEY);
                        const bool tntHaveSub = (tntSubYHi > tntSubYLo) &&
                        (!tntSubTDX.empty() || !tntSubNTDX.empty());
                        const double tntSubPad = tntHaveSub ? std::max(0.35 * (tntSubYHi - tntSubYLo), 0.5) : 1.0;
                        
                        double tntSigYLo = 1e30, tntSigYHi = -1e30;
                        auto TntUpdateSigRange = [&](const vector<double>& y, const vector<double>& ey) {
                            for (std::size_t i = 0; i < y.size(); ++i) {
                                tntSigYLo = std::min(tntSigYLo, y[i] - ey[i]);
                                tntSigYHi = std::max(tntSigYHi, y[i] + ey[i]); } };
                        TntUpdateSigRange(tntSubTDSigY, tntSubTDSigEY);
                        TntUpdateSigRange(tntSubNTDSigY, tntSubNTDSigEY);
                        const double tntSigPad = tntHaveSub ? std::max(0.35 * (tntSigYHi - tntSigYLo), 0.5) : 1.0;
                        
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
                                         "c_tightNonTight", 900, tntHaveSub ? 1000 : 700
                                         );
                            cTNT.cd();
                            
                            // Upper pad
                            const double tntPadLoEdge = tntHaveSub ? 0.36 : 0.0;
                            TPad* padUpTNT = new TPad("padUpTNT", "padUpTNT", 0.0, tntPadLoEdge, 1.0, 1.0);
                            padUpTNT->SetBottomMargin(tntHaveSub ? 0.10 : 0.12);
                            padUpTNT->SetLeftMargin(0.14);
                            padUpTNT->SetRightMargin(0.04);
                            padUpTNT->SetTopMargin(0.08);
                            padUpTNT->Draw();
                            padUpTNT->cd();
                            
                            hTMc->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                            hTMc->GetYaxis()->SetTitle("Normalized to unit area");
                            hTMc->GetXaxis()->SetTitleSize(0.045);
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
                            
                            TF1* fTDataGauss  = DrawGaussFitCurve(hTData,  kBlack);
                            TF1* fNTDataGauss = DrawGaussFitCurve(hNTData, kRed+1);
                            
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
                            
                            // Middle pad: Gaussian mean vs pT
                            if (tntHaveSub)
                            {
                                cTNT.cd();
                                TPad* padMidTNT = new TPad("padMidTNT", "padMidTNT", 0.0, 0.21, 1.0, 0.36);
                                padMidTNT->SetTopMargin(0.02);
                                padMidTNT->SetBottomMargin(0.00);
                                padMidTNT->SetLeftMargin(0.14);
                                padMidTNT->SetRightMargin(0.04);
                                padMidTNT->Draw();
                                padMidTNT->cd();
                                
                                TH1F* hFrSubTNT = new TH1F(
                                                           TString::Format("hFrSub_tnt_%s_%s_%s_%s",
                                                                           mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                                           "", 100, 10.0, kPtEdges.back());
                                hFrSubTNT->SetDirectory(nullptr); hFrSubTNT->SetStats(0);
                                hFrSubTNT->SetMinimum(tntSubYLo - tntSubPad);
                                hFrSubTNT->SetMaximum(tntSubYHi + tntSubPad);
                                hFrSubTNT->GetYaxis()->SetTitle("#mu^{Gauss}[GeV]");
                                hFrSubTNT->GetYaxis()->SetTitleSize(0.19);
                                hFrSubTNT->GetYaxis()->SetTitleOffset(0.28);
                                hFrSubTNT->GetYaxis()->SetLabelSize(0.14);
                                hFrSubTNT->GetYaxis()->SetNdivisions(505);
                                hFrSubTNT->GetXaxis()->SetTitle("");
                                hFrSubTNT->GetXaxis()->SetTitleSize(0.0);
                                hFrSubTNT->GetXaxis()->SetLabelSize(0.0);
                                hFrSubTNT->GetXaxis()->SetTickLength(0.0);
                                hFrSubTNT->Draw();
                                
                                auto MakeSubG = [](const vector<double>& x, const vector<double>& y, const vector<double>& ey,
                                                   int marker, int color) -> TGraphErrors* {
                                    if (x.empty()) return nullptr;
                                    vector<double> ex(x.size(), 0.0);
                                    TGraphErrors* g = new TGraphErrors((int)x.size(), &x[0], &y[0], &ex[0], &ey[0]);
                                    g->SetMarkerStyle(marker); g->SetMarkerSize(1.0);
                                    g->SetMarkerColor(color);  g->SetLineColor(color);
                                    g->SetLineWidth(2); g->Draw("PE1 SAME"); return g; };
                                
                                TGraphErrors* gTD  = MakeSubG(tntSubTDX,  tntSubTDY,  tntSubTDEY,  20, kBlack);
                                TGraphErrors* gNTD = MakeSubG(tntSubNTDX, tntSubNTDY, tntSubNTDEY, 24, kRed+1);
                                
                                // Bottom pad: Gaussian sigma vs pT
                                cTNT.cd();
                                TPad* padBotTNT = new TPad("padBotTNT", "padBotTNT", 0.0, 0.0, 1.0, 0.21);
                                padBotTNT->SetTopMargin(0.00);
                                padBotTNT->SetBottomMargin(0.30);
                                padBotTNT->SetLeftMargin(0.14);
                                padBotTNT->SetRightMargin(0.04);
                                padBotTNT->Draw();
                                padBotTNT->cd();
                                
                                TH1F* hFrSigTNT = new TH1F(
                                                           TString::Format("hFrSig_tnt_%s_%s_%s_%s",
                                                                           mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                                           "", 100, 10.0, kPtEdges.back());
                                hFrSigTNT->SetDirectory(nullptr); hFrSigTNT->SetStats(0);
                                hFrSigTNT->SetMinimum(tntSigYLo - tntSigPad);
                                hFrSigTNT->SetMaximum(tntSigYHi + tntSigPad);
                                hFrSigTNT->GetYaxis()->SetTitle("#sigma^{Gauss}[GeV]");
                                hFrSigTNT->GetYaxis()->SetTitleSize(0.14);
                                hFrSigTNT->GetYaxis()->SetTitleOffset(0.32);
                                hFrSigTNT->GetYaxis()->SetLabelSize(0.10);
                                hFrSigTNT->GetYaxis()->SetNdivisions(505);
                                hFrSigTNT->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                                hFrSigTNT->GetXaxis()->SetTitleSize(0.12);
                                hFrSigTNT->GetXaxis()->SetTitleOffset(0.95);
                                hFrSigTNT->GetXaxis()->SetLabelSize(0.10);
                                hFrSigTNT->Draw();
                                
                                TGraphErrors* gTDSig  = MakeSubG(tntSubTDX,  tntSubTDSigY,  tntSubTDSigEY,  20, kBlack);
                                TGraphErrors* gNTDSig = MakeSubG(tntSubNTDX, tntSubNTDSigY, tntSubNTDSigEY, 24, kRed+1);
                                
                                cTNT.Modified();
                                cTNT.Update();
                                
                                SaveCanvas(cTNT, JoinPath(ptDir, "Eiso_tightNonTight_overlay.png"));
                                
                                if (gTDSig)  delete gTDSig;
                                if (gNTDSig) delete gNTDSig;
                                delete hFrSigTNT;
                                if (gTD)  delete gTD;
                                if (gNTD) delete gNTD;
                                delete hFrSubTNT;
                            }
                            
                            if (!tntHaveSub)
                                SaveCanvas(cTNT, JoinPath(ptDir, "Eiso_tightNonTight_overlay.png"));
                            
                            if (fTDataGauss)  delete fTDataGauss;
                            if (fNTDataGauss) delete fNTDataGauss;
                            
                            // --- save 4 individual unnormalized distributions ---s
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
                        
                        {
                            auto UpdateSummaryRange = [&](const vector<double>& y, const vector<double>& ey, double& yLo, double& yHi)
                            {
                                for (std::size_t i = 0; i < y.size(); ++i)
                                {
                                    yLo = std::min(yLo, y[i] - ey[i]);
                                    yHi = std::max(yHi, y[i] + ey[i]);
                                }
                            };
                            
                            auto MakeSummaryG = [](const vector<double>& x, const vector<double>& y, const vector<double>& ey,
                                                   int marker, int color) -> TGraphErrors*
                            {
                                if (x.empty()) return nullptr;
                                vector<double> ex(x.size(), 0.0);
                                TGraphErrors* g = new TGraphErrors((int)x.size(), &x[0], &y[0], &ex[0], &ey[0]);
                                g->SetMarkerStyle(marker);
                                g->SetMarkerSize(1.1);
                                g->SetMarkerColor(color);
                                g->SetLineColor(color);
                                g->SetLineWidth(2);
                                return g;
                            };
                            
                            const double coneRValTNT = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                            
                            double tntMeanYLo = 1e30, tntMeanYHi = -1e30;
                            UpdateSummaryRange(tntSubTDY,   tntSubTDEY,   tntMeanYLo, tntMeanYHi);
                            UpdateSummaryRange(tntSubNTDY,  tntSubNTDEY,  tntMeanYLo, tntMeanYHi);
                            UpdateSummaryRange(tntSubTMY,   tntSubTMEY,   tntMeanYLo, tntMeanYHi);
                            UpdateSummaryRange(tntSubNTMY,  tntSubNTMEY,  tntMeanYLo, tntMeanYHi);
                            const bool haveTntMeanSummary =
                            (!tntSubTDX.empty() || !tntSubNTDX.empty() || !tntSubTMX.empty() || !tntSubNTMX.empty());
                            
                            if (haveTntMeanSummary)
                            {
                                const double tntMeanPad = (tntMeanYHi > tntMeanYLo) ? std::max(0.25 * (tntMeanYHi - tntMeanYLo), 0.25) : 0.25;
                                
                                TCanvas cMeanTNT(
                                                 TString::Format("cMeanTNT_%s_%s_%s", mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str()).Data(),
                                                 "cMeanTNT", 900, 700);
                                ApplyCanvasMargins1D(cMeanTNT);
                                cMeanTNT.cd();
                                
                                TH1F hFrMeanTNT(
                                               TString::Format("hFrMeanTNT_%s_%s_%s", mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str()).Data(),
                                               "", 100, 10.0, kPtEdges.back());
                                hFrMeanTNT.SetDirectory(nullptr);
                                hFrMeanTNT.SetStats(0);
                                hFrMeanTNT.SetMinimum(std::max(0.0, tntMeanYLo - tntMeanPad));
                                hFrMeanTNT.SetMaximum(tntMeanYHi + tntMeanPad);
                                hFrMeanTNT.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                                hFrMeanTNT.GetYaxis()->SetTitle("#mu^{Gauss} [GeV]");
                                hFrMeanTNT.GetXaxis()->SetTitleSize(0.055);
                                hFrMeanTNT.GetYaxis()->SetTitleSize(0.055);
                                hFrMeanTNT.GetXaxis()->SetLabelSize(0.045);
                                hFrMeanTNT.GetYaxis()->SetLabelSize(0.045);
                                hFrMeanTNT.GetYaxis()->SetTitleOffset(1.15);
                                hFrMeanTNT.Draw();
                                
                                TGraphErrors* gMeanTD   = MakeSummaryG(tntSubTDX,   tntSubTDY,   tntSubTDEY,   20, kBlack);
                                TGraphErrors* gMeanNTD  = MakeSummaryG(tntSubNTDX,  tntSubNTDY,  tntSubNTDEY,  24, kRed + 1);
                                TGraphErrors* gMeanTM   = MakeSummaryG(tntSubTMX,   tntSubTMY,   tntSubTMEY,   25, kBlue + 1);
                                TGraphErrors* gMeanNTM  = MakeSummaryG(tntSubNTMX,  tntSubNTMY,  tntSubNTMEY,  26, kGreen + 2);
                                
                                if (gMeanTD)  gMeanTD->Draw("PE1 SAME");
                                if (gMeanNTD) gMeanNTD->Draw("PE1 SAME");
                                if (gMeanTM)  gMeanTM->Draw("PE1 SAME");
                                if (gMeanNTM) gMeanNTM->Draw("PE1 SAME");
                                
                                TLegend legMeanTNT(0.48, 0.66, 0.92, 0.88);
                                legMeanTNT.SetBorderSize(0);
                                legMeanTNT.SetFillStyle(0);
                                legMeanTNT.SetTextFont(42);
                                legMeanTNT.SetTextSize(0.032);
                                legMeanTNT.SetNColumns(2);
                                if (gMeanTD)  legMeanTNT.AddEntry(gMeanTD,  "tight data", "ep");
                                if (gMeanNTD) legMeanTNT.AddEntry(gMeanNTD, "nontight data", "ep");
                                if (gMeanTM)  legMeanTNT.AddEntry(gMeanTM,  tightMcLegend.c_str(), "ep");
                                if (gMeanNTM) legMeanTNT.AddEntry(gMeanNTM, nonTightMcLegend.c_str(), "ep");
                                legMeanTNT.Draw();
                                
                                TLatex tMeanTitle;
                                tMeanTitle.SetNDC(true);
                                tMeanTitle.SetTextFont(42);
                                tMeanTitle.SetTextAlign(23);
                                tMeanTitle.SetTextSize(0.038);
                                tMeanTitle.DrawLatex(0.50, 0.97,
                                                     TString::Format("#mu^{Gauss}: tight/nontight Au+Au data vs %s embedded MC",
                                                                     mcTitleTag.c_str()).Data());
                                
                                TLatex tMeanInfo;
                                tMeanInfo.SetNDC(true);
                                tMeanInfo.SetTextFont(42);
                                tMeanInfo.SetTextAlign(13);
                                tMeanInfo.SetTextSize(0.034);
                                tMeanInfo.DrawLatex(0.18, 0.88, trigDisplayLabel.c_str());
                                tMeanInfo.DrawLatex(0.18, 0.84, TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                                tMeanInfo.DrawLatex(0.18, 0.80, TString::Format("UE: %s", dataH.label.c_str()).Data());
                                tMeanInfo.DrawLatex(0.18, 0.76, TString::Format("#DeltaR_{cone} < %.1f", coneRValTNT).Data());
                                
                                TLatex tMeanSph;
                                tMeanSph.SetNDC(true);
                                tMeanSph.SetTextFont(42);
                                tMeanSph.SetTextAlign(33);
                                tMeanSph.SetTextSize(0.042);
                                tMeanSph.DrawLatex(0.92, 0.18, "#bf{sPHENIX} #it{Internal}");
                                tMeanSph.SetTextSize(0.034);
                                tMeanSph.DrawLatex(0.92, 0.13, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                                
                                SaveCanvas(cMeanTNT, JoinPath(meanSummaryDir,
                                                              TString::Format("tightNonTight_meanOverlay_%s.png", cb.folder.c_str()).Data()));
                                
                                if (mcOverlayFolder == "photonJetOverlays" && (gMeanTD || gMeanNTD))
                                {
                                    TCanvas cMeanTNTDataOnly(
                                                             TString::Format("cMeanTNTDataOnly_%s_%s_%s", mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str()).Data(),
                                                             "cMeanTNTDataOnly", 900, 700);
                                    ApplyCanvasMargins1D(cMeanTNTDataOnly);
                                    cMeanTNTDataOnly.cd();
                                    
                                    TH1F hFrMeanTNTDataOnly(
                                                           TString::Format("hFrMeanTNTDataOnly_%s_%s_%s", mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str()).Data(),
                                                           "", 100, 10.0, kPtEdges.back());
                                    hFrMeanTNTDataOnly.SetDirectory(nullptr);
                                    hFrMeanTNTDataOnly.SetStats(0);
                                    hFrMeanTNTDataOnly.SetMinimum(std::max(0.0, tntMeanYLo - tntMeanPad));
                                    hFrMeanTNTDataOnly.SetMaximum(tntMeanYHi + tntMeanPad);
                                    hFrMeanTNTDataOnly.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                                    hFrMeanTNTDataOnly.GetYaxis()->SetTitle("#mu^{Gauss} [GeV]");
                                    hFrMeanTNTDataOnly.GetXaxis()->SetTitleSize(0.055);
                                    hFrMeanTNTDataOnly.GetYaxis()->SetTitleSize(0.055);
                                    hFrMeanTNTDataOnly.GetXaxis()->SetLabelSize(0.045);
                                    hFrMeanTNTDataOnly.GetYaxis()->SetLabelSize(0.045);
                                    hFrMeanTNTDataOnly.GetYaxis()->SetTitleOffset(1.15);
                                    hFrMeanTNTDataOnly.Draw();
                                    
                                    if (gMeanTD)  gMeanTD->Draw("PE1 SAME");
                                    if (gMeanNTD) gMeanNTD->Draw("PE1 SAME");
                                    
                                    TLegend legMeanTNTDataOnly(0.56, 0.74, 0.92, 0.88);
                                    legMeanTNTDataOnly.SetBorderSize(0);
                                    legMeanTNTDataOnly.SetFillStyle(0);
                                    legMeanTNTDataOnly.SetTextFont(42);
                                    legMeanTNTDataOnly.SetTextSize(0.032);
                                    if (gMeanTD)  legMeanTNTDataOnly.AddEntry(gMeanTD,  "tight data", "ep");
                                    if (gMeanNTD) legMeanTNTDataOnly.AddEntry(gMeanNTD, "nontight data", "ep");
                                    legMeanTNTDataOnly.Draw();
                                    
                                    TLatex tMeanTitleDataOnly;
                                    tMeanTitleDataOnly.SetNDC(true);
                                    tMeanTitleDataOnly.SetTextFont(42);
                                    tMeanTitleDataOnly.SetTextAlign(23);
                                    tMeanTitleDataOnly.SetTextSize(0.038);
                                    tMeanTitleDataOnly.DrawLatex(0.50, 0.97,
                                                                 "#mu^{Gauss}: tight/nontight Au+Au data only");
                                    
                                    TLatex tMeanInfoDataOnly;
                                    tMeanInfoDataOnly.SetNDC(true);
                                    tMeanInfoDataOnly.SetTextFont(42);
                                    tMeanInfoDataOnly.SetTextAlign(13);
                                    tMeanInfoDataOnly.SetTextSize(0.034);
                                    tMeanInfoDataOnly.DrawLatex(0.18, 0.88, trigDisplayLabel.c_str());
                                    tMeanInfoDataOnly.DrawLatex(0.18, 0.84, TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                                    tMeanInfoDataOnly.DrawLatex(0.18, 0.80, TString::Format("UE: %s", dataH.label.c_str()).Data());
                                    tMeanInfoDataOnly.DrawLatex(0.18, 0.76, TString::Format("#DeltaR_{cone} < %.1f", coneRValTNT).Data());
                                    
                                    TLatex tMeanSphDataOnly;
                                    tMeanSphDataOnly.SetNDC(true);
                                    tMeanSphDataOnly.SetTextFont(42);
                                    tMeanSphDataOnly.SetTextAlign(33);
                                    tMeanSphDataOnly.SetTextSize(0.042);
                                    tMeanSphDataOnly.DrawLatex(0.92, 0.18, "#bf{sPHENIX} #it{Internal}");
                                    tMeanSphDataOnly.SetTextSize(0.034);
                                    tMeanSphDataOnly.DrawLatex(0.92, 0.13, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                                    
                                    SaveCanvas(cMeanTNTDataOnly, JoinPath(meanSummaryDir,
                                                                          TString::Format("tightNonTight_meanOverlay_dataOnly_%s.png", cb.folder.c_str()).Data()));
                                }
                                
                                if (gMeanTD) delete gMeanTD;
                                if (gMeanNTD) delete gMeanNTD;
                                if (gMeanTM) delete gMeanTM;
                                if (gMeanNTM) delete gMeanNTM;
                            }
                            
                            double tntSigmaYLo = 1e30, tntSigmaYHi = -1e30;
                            UpdateSummaryRange(tntSubTDSigY,   tntSubTDSigEY,   tntSigmaYLo, tntSigmaYHi);
                            UpdateSummaryRange(tntSubNTDSigY,  tntSubNTDSigEY,  tntSigmaYLo, tntSigmaYHi);
                            UpdateSummaryRange(tntSubTMSigY,   tntSubTMSigEY,   tntSigmaYLo, tntSigmaYHi);
                            UpdateSummaryRange(tntSubNTMSigY,  tntSubNTMSigEY,  tntSigmaYLo, tntSigmaYHi);
                            const bool haveTntSigmaSummary =
                            (!tntSubTDX.empty() || !tntSubNTDX.empty() || !tntSubTMX.empty() || !tntSubNTMX.empty());
                            
                            if (haveTntSigmaSummary)
                            {
                                const double tntSigmaPad = (tntSigmaYHi > tntSigmaYLo) ? std::max(0.25 * (tntSigmaYHi - tntSigmaYLo), 0.25) : 0.25;
                                
                                TCanvas cSigmaTNT(
                                                  TString::Format("cSigmaTNT_%s_%s_%s", mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str()).Data(),
                                                  "cSigmaTNT", 900, 700);
                                ApplyCanvasMargins1D(cSigmaTNT);
                                cSigmaTNT.cd();
                                
                                TH1F hFrSigmaTNT(
                                                TString::Format("hFrSigmaTNT_%s_%s_%s", mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str()).Data(),
                                                "", 100, 10.0, kPtEdges.back());
                                hFrSigmaTNT.SetDirectory(nullptr);
                                hFrSigmaTNT.SetStats(0);
                                hFrSigmaTNT.SetMinimum(std::max(0.0, tntSigmaYLo - tntSigmaPad));
                                hFrSigmaTNT.SetMaximum(tntSigmaYHi + tntSigmaPad);
                                hFrSigmaTNT.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                                hFrSigmaTNT.GetYaxis()->SetTitle("#sigma^{Gauss} [GeV]");
                                hFrSigmaTNT.GetXaxis()->SetTitleSize(0.055);
                                hFrSigmaTNT.GetYaxis()->SetTitleSize(0.055);
                                hFrSigmaTNT.GetXaxis()->SetLabelSize(0.045);
                                hFrSigmaTNT.GetYaxis()->SetLabelSize(0.045);
                                hFrSigmaTNT.GetYaxis()->SetTitleOffset(1.15);
                                hFrSigmaTNT.Draw();
                                
                                TGraphErrors* gSigmaTD   = MakeSummaryG(tntSubTDX,   tntSubTDSigY,   tntSubTDSigEY,   20, kBlack);
                                TGraphErrors* gSigmaNTD  = MakeSummaryG(tntSubNTDX,  tntSubNTDSigY,  tntSubNTDSigEY,  24, kRed + 1);
                                TGraphErrors* gSigmaTM   = MakeSummaryG(tntSubTMX,   tntSubTMSigY,   tntSubTMSigEY,   25, kBlue + 1);
                                TGraphErrors* gSigmaNTM  = MakeSummaryG(tntSubNTMX,  tntSubNTMSigY,  tntSubNTMSigEY,  26, kGreen + 2);
                                
                                if (gSigmaTD)  gSigmaTD->Draw("PE1 SAME");
                                if (gSigmaNTD) gSigmaNTD->Draw("PE1 SAME");
                                if (gSigmaTM)  gSigmaTM->Draw("PE1 SAME");
                                if (gSigmaNTM) gSigmaNTM->Draw("PE1 SAME");
                                
                                TLegend legSigmaTNT(0.48, 0.66, 0.92, 0.88);
                                legSigmaTNT.SetBorderSize(0);
                                legSigmaTNT.SetFillStyle(0);
                                legSigmaTNT.SetTextFont(42);
                                legSigmaTNT.SetTextSize(0.032);
                                legSigmaTNT.SetNColumns(2);
                                if (gSigmaTD)  legSigmaTNT.AddEntry(gSigmaTD,  "tight data", "ep");
                                if (gSigmaNTD) legSigmaTNT.AddEntry(gSigmaNTD, "nontight data", "ep");
                                if (gSigmaTM)  legSigmaTNT.AddEntry(gSigmaTM,  tightMcLegend.c_str(), "ep");
                                if (gSigmaNTM) legSigmaTNT.AddEntry(gSigmaNTM, nonTightMcLegend.c_str(), "ep");
                                legSigmaTNT.Draw();
                                
                                TLatex tSigmaTitle;
                                tSigmaTitle.SetNDC(true);
                                tSigmaTitle.SetTextFont(42);
                                tSigmaTitle.SetTextAlign(23);
                                tSigmaTitle.SetTextSize(0.038);
                                tSigmaTitle.DrawLatex(0.50, 0.97,
                                                      TString::Format("#sigma^{Gauss}: tight/nontight Au+Au data vs %s embedded MC",
                                                                      mcTitleTag.c_str()).Data());
                                
                                TLatex tSigmaInfo;
                                tSigmaInfo.SetNDC(true);
                                tSigmaInfo.SetTextFont(42);
                                tSigmaInfo.SetTextAlign(13);
                                tSigmaInfo.SetTextSize(0.034);
                                tSigmaInfo.DrawLatex(0.18, 0.88, trigDisplayLabel.c_str());
                                tSigmaInfo.DrawLatex(0.18, 0.84, TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                                tSigmaInfo.DrawLatex(0.18, 0.80, TString::Format("UE: %s", dataH.label.c_str()).Data());
                                tSigmaInfo.DrawLatex(0.18, 0.76, TString::Format("#DeltaR_{cone} < %.1f", coneRValTNT).Data());
                                
                                TLatex tSigmaSph;
                                tSigmaSph.SetNDC(true);
                                tSigmaSph.SetTextFont(42);
                                tSigmaSph.SetTextAlign(33);
                                tSigmaSph.SetTextSize(0.042);
                                tSigmaSph.DrawLatex(0.92, 0.18, "#bf{sPHENIX} #it{Internal}");
                                tSigmaSph.SetTextSize(0.034);
                                tSigmaSph.DrawLatex(0.92, 0.13, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                                
                                SaveCanvas(cSigmaTNT, JoinPath(sigmaSummaryDir,
                                                               TString::Format("tightNonTight_sigmaOverlay_%s.png", cb.folder.c_str()).Data()));
                                
                                if (mcOverlayFolder == "photonJetOverlays" && (gSigmaTD || gSigmaNTD))
                                {
                                    TCanvas cSigmaTNTDataOnly(
                                                              TString::Format("cSigmaTNTDataOnly_%s_%s_%s", mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str()).Data(),
                                                              "cSigmaTNTDataOnly", 900, 700);
                                    ApplyCanvasMargins1D(cSigmaTNTDataOnly);
                                    cSigmaTNTDataOnly.cd();
                                    
                                    TH1F hFrSigmaTNTDataOnly(
                                                             TString::Format("hFrSigmaTNTDataOnly_%s_%s_%s", mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str()).Data(),
                                                             "", 100, 10.0, kPtEdges.back());
                                    hFrSigmaTNTDataOnly.SetDirectory(nullptr);
                                    hFrSigmaTNTDataOnly.SetStats(0);
                                    hFrSigmaTNTDataOnly.SetMinimum(std::max(0.0, tntSigmaYLo - tntSigmaPad));
                                    hFrSigmaTNTDataOnly.SetMaximum(tntSigmaYHi + tntSigmaPad);
                                    hFrSigmaTNTDataOnly.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                                    hFrSigmaTNTDataOnly.GetYaxis()->SetTitle("#sigma^{Gauss} [GeV]");
                                    hFrSigmaTNTDataOnly.GetXaxis()->SetTitleSize(0.055);
                                    hFrSigmaTNTDataOnly.GetYaxis()->SetTitleSize(0.055);
                                    hFrSigmaTNTDataOnly.GetXaxis()->SetLabelSize(0.045);
                                    hFrSigmaTNTDataOnly.GetYaxis()->SetLabelSize(0.045);
                                    hFrSigmaTNTDataOnly.GetYaxis()->SetTitleOffset(1.15);
                                    hFrSigmaTNTDataOnly.Draw();
                                    
                                    if (gSigmaTD)  gSigmaTD->Draw("PE1 SAME");
                                    if (gSigmaNTD) gSigmaNTD->Draw("PE1 SAME");
                                    
                                    TLegend legSigmaTNTDataOnly(0.56, 0.74, 0.92, 0.88);
                                    legSigmaTNTDataOnly.SetBorderSize(0);
                                    legSigmaTNTDataOnly.SetFillStyle(0);
                                    legSigmaTNTDataOnly.SetTextFont(42);
                                    legSigmaTNTDataOnly.SetTextSize(0.032);
                                    if (gSigmaTD)  legSigmaTNTDataOnly.AddEntry(gSigmaTD,  "tight data", "ep");
                                    if (gSigmaNTD) legSigmaTNTDataOnly.AddEntry(gSigmaNTD, "nontight data", "ep");
                                    legSigmaTNTDataOnly.Draw();
                                    
                                    TLatex tSigmaTitleDataOnly;
                                    tSigmaTitleDataOnly.SetNDC(true);
                                    tSigmaTitleDataOnly.SetTextFont(42);
                                    tSigmaTitleDataOnly.SetTextAlign(23);
                                    tSigmaTitleDataOnly.SetTextSize(0.038);
                                    tSigmaTitleDataOnly.DrawLatex(0.50, 0.97,
                                                                  "#sigma^{Gauss}: tight/nontight Au+Au data only");
                                    
                                    TLatex tSigmaInfoDataOnly;
                                    tSigmaInfoDataOnly.SetNDC(true);
                                    tSigmaInfoDataOnly.SetTextFont(42);
                                    tSigmaInfoDataOnly.SetTextAlign(13);
                                    tSigmaInfoDataOnly.SetTextSize(0.034);
                                    tSigmaInfoDataOnly.DrawLatex(0.18, 0.88, trigDisplayLabel.c_str());
                                    tSigmaInfoDataOnly.DrawLatex(0.18, 0.84, TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                                    tSigmaInfoDataOnly.DrawLatex(0.18, 0.80, TString::Format("UE: %s", dataH.label.c_str()).Data());
                                    tSigmaInfoDataOnly.DrawLatex(0.18, 0.76, TString::Format("#DeltaR_{cone} < %.1f", coneRValTNT).Data());
                                    
                                    TLatex tSigmaSphDataOnly;
                                    tSigmaSphDataOnly.SetNDC(true);
                                    tSigmaSphDataOnly.SetTextFont(42);
                                    tSigmaSphDataOnly.SetTextAlign(33);
                                    tSigmaSphDataOnly.SetTextSize(0.042);
                                    tSigmaSphDataOnly.DrawLatex(0.92, 0.18, "#bf{sPHENIX} #it{Internal}");
                                    tSigmaSphDataOnly.SetTextSize(0.034);
                                    tSigmaSphDataOnly.DrawLatex(0.92, 0.13, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                                    
                                    SaveCanvas(cSigmaTNTDataOnly, JoinPath(sigmaSummaryDir,
                                                                           TString::Format("tightNonTight_sigmaOverlay_dataOnly_%s.png", cb.folder.c_str()).Data()));
                                }
                                
                                if (gSigmaTD) delete gSigmaTD;
                                if (gSigmaNTD) delete gSigmaNTD;
                                if (gSigmaTM) delete gSigmaTM;
                                if (gSigmaNTM) delete gSigmaNTM;
                            }
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
    if (fPPSimMerged)
    {
        fPPSimMerged->Close();
        delete fPPSimMerged;
        fPPSimMerged = nullptr;
    }
    if (fPP)
    {
        fPP->Close();
        delete fPP;
        fPP = nullptr;
    }
}
