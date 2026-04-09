// SS_QA sits alongside isoQA under the same trigger output dir
const string ssQADir  = JoinPath(trigOutBase, "SS_QA");
EnsureDir(ssQADir);

struct SSVarDef { string var; string label; };
const vector<SSVarDef> ssVars = {
  {"weta",   "w_{#eta}"},
  {"wphi",   "w_{#phi}"},
  {"e11e33", "E_{11}/E_{33}"},
  {"et1",    "et1"},
  {"e32e35", "E_{32}/E_{35}"}
};

auto TagLabel = [&](const string& tag) -> string
{
  if (tag == "inclusive") return "Before preselection";
  if (tag == "pre") return "Preselection";
  if (tag == "tight") return "Tight";
  if (tag == "nonTight") return "Non-tight";
  return tag;
};

auto TagFolder = [&](const string& tag) -> string
{
  if (tag == "inclusive") return "beforePreselection";
  if (tag == "pre") return "preselection";
  if (tag == "tight") return "tight";
  if (tag == "nonTight") return "nonTight";
  return tag;
};

const vector<string> ppg12Tags = {"inclusive", "pre", "tight", "nonTight"};

const int ssQA_rebinFactor = 3;
std::set<TH1*> ssQA_alreadyRebinned;

auto GetTH1FromTopDir = [&](TDirectory* topDir, const string& hname) -> TH1*
{
  if (!topDir) return nullptr;
  TObject* obj = topDir->Get(hname.c_str());
  if (!obj) return nullptr;
  TH1* h = dynamic_cast<TH1*>(obj);
  if (h && ssQA_rebinFactor > 1 && ssQA_alreadyRebinned.find(h) == ssQA_alreadyRebinned.end())
  {
    h->Rebin(ssQA_rebinFactor);
    ssQA_alreadyRebinned.insert(h);
  }
  return h;
};

auto StyleOverlayHist = [&](TH1* h, int color, int mstyle) -> void
{
  if (!h) return;
  h->SetLineWidth(2);
  h->SetLineColor(color);
  h->SetMarkerStyle(mstyle);
  h->SetMarkerSize(0.95);
  h->SetMarkerColor(color);
};

auto DrawMissingPad = [&](const string& titleLine) -> void
{
  TLatex t;
  t.SetNDC(true);
  t.SetTextFont(42);
  t.SetTextAlign(22);
  t.SetTextSize(0.080);
  t.DrawLatex(0.50, 0.55, "MISSING");

  t.SetTextSize(0.050);
  t.DrawLatex(0.50, 0.42, titleLine.c_str());
};

auto CloneNormalizeStyle = [&](TH1* hIn,
                               const string& newName,
                               int color,
                               int mstyle) -> TH1*
{
  if (!hIn) return nullptr;

  TH1* h = CloneTH1(hIn, newName);
  if (!h) return nullptr;

  EnsureSumw2(h);

  const double I = h->Integral(0, h->GetNbinsX() + 1);
  if (I > 0.0) h->Scale(1.0 / I);

  StyleOverlayHist(h, color, mstyle);
  h->SetTitle("");

  return h;
};

auto PassDefValidForVar = [&](const string& var, const string& passDef) -> bool
{
  if (passDef == "pre")
  {
    return (var == "weta" || var == "e11e33" || var == "et1" || var == "e32e35");
  }
  if (passDef == "tight")
  {
    return (var == "weta" || var == "wphi" || var == "e11e33" || var == "et1" || var == "e32e35");
  }
  return false;
};

auto PassFractionForHist = [&](TH1* h,
                               const string& var,
                               const string& passDef,
                               double ptCenter,
                               double& frac,
                               double& err) -> bool
{
  frac = 0.0;
  err  = 0.0;

  if (!h) return false;

  const int nb = h->GetNbinsX();
  const double den = h->Integral(0, nb + 1);
  if (!(den > 0.0)) return false;

  double num = 0.0;

  if (var == "weta")
  {
    if (passDef == "pre")
    {
      const double cut = 0.6;
      const int bHi = h->GetXaxis()->FindBin(cut - 1e-6);
      num = h->Integral(1, bHi);
    }
    else if (passDef == "tight")
    {
      const double cutLo = 0.0;
      const double cutHi = 0.15 + 0.006 * ptCenter;
      const int bLo = h->GetXaxis()->FindBin(cutLo + 1e-6);
      const int bHi = h->GetXaxis()->FindBin(cutHi - 1e-6);
      num = h->Integral(bLo, bHi);
    }
    else return false;
  }
  else if (var == "wphi")
  {
    if (passDef == "tight")
    {
      const double cutLo = 0.0;
      const double cutHi = 0.15 + 0.006 * ptCenter;
      const int bLo = h->GetXaxis()->FindBin(cutLo + 1e-6);
      const int bHi = h->GetXaxis()->FindBin(cutHi - 1e-6);
      num = h->Integral(bLo, bHi);
    }
    else return false;
  }
  else if (var == "e11e33")
  {
    if (passDef == "pre")
    {
      const double cutHi = 0.98;
      const int bHi = h->GetXaxis()->FindBin(cutHi - 1e-6);
      num = h->Integral(1, bHi);
    }
    else if (passDef == "tight")
    {
      const double cutLo = 0.4;
      const double cutHi = 0.98;
      const int bLo = h->GetXaxis()->FindBin(cutLo + 1e-6);
      const int bHi = h->GetXaxis()->FindBin(cutHi - 1e-6);
      num = h->Integral(bLo, bHi);
    }
    else return false;
  }
  else if (var == "et1")
  {
    if (passDef == "pre")
    {
      const double cutLo = 0.6;
      const double cutHi = 1.0;
      const int bLo = h->GetXaxis()->FindBin(cutLo + 1e-6);
      const int bHi = h->GetXaxis()->FindBin(cutHi - 1e-6);
      num = h->Integral(bLo, bHi);
    }
    else if (passDef == "tight")
    {
      const double cutLo = 0.9;
      const double cutHi = 1.0;
      const int bLo = h->GetXaxis()->FindBin(cutLo + 1e-6);
      const int bHi = h->GetXaxis()->FindBin(cutHi - 1e-6);
      num = h->Integral(bLo, bHi);
    }
    else return false;
  }
  else if (var == "e32e35")
  {
    if (passDef == "pre")
    {
      const double cutLo = 0.8;
      const double cutHi = 1.0;
      const int bLo = h->GetXaxis()->FindBin(cutLo + 1e-6);
      const int bHi = h->GetXaxis()->FindBin(cutHi - 1e-6);
      num = h->Integral(bLo, bHi);
    }
    else if (passDef == "tight")
    {
      const double cutLo = 0.92;
      const double cutHi = 1.0;
      const int bLo = h->GetXaxis()->FindBin(cutLo + 1e-6);
      const int bHi = h->GetXaxis()->FindBin(cutHi - 1e-6);
      num = h->Integral(bLo, bHi);
    }
    else return false;
  }
  else
  {
    return false;
  }

  frac = (den > 0.0) ? (num / den) : 0.0;
  if (frac < 0.0) frac = 0.0;
  if (frac > 1.0) frac = 1.0;

  err = std::sqrt(frac * std::max(0.0, 1.0 - frac) / den);
  return true;
};

// UE variant indices: 0 = noSub, 2 = variantA, 3 = variantB (skip 1 = baseVariant)
const vector<std::size_t> ssUEIndices = {std::size_t(0), std::size_t(2), std::size_t(3)};
const int ssColors[4] = {kBlack, kBlue + 1, kOrange + 7, kGreen + 2};
const int ssMarkers[4] = {20, 20, 20, 20};

// Force the SS template SIM source to the merged photonJet5+10+20 file,
// independent of the current SIM toggle state.
TFile* fSimSS = nullptr;
TDirectory* simTopSS = nullptr;
{
  const string simMergedSS =
    MergedSimPath("photonJet5and10and20merged_SIM", "RecoilJets_photonjet5plus10plus20_MERGED.root");

  bool haveMergedSS = false;
  if (!simMergedSS.empty() && !gSystem->AccessPathName(simMergedSS.c_str()))
  {
    haveMergedSS = true;
  }
  else
  {
    const bool okMergeSS = BuildMergedSIMFile_PhotonSlices(
      {InputSim("photonjet5"), InputSim("photonjet10"), InputSim("photonjet20")},
      {kSigmaPhoton5_pb, kSigmaPhoton10_pb, kSigmaPhoton20_pb},
      simMergedSS,
      kDirSIM,
      {"photonJet5", "photonJet10", "photonJet20"}
    );
    haveMergedSS = okMergeSS && !gSystem->AccessPathName(simMergedSS.c_str());
  }

  if (haveMergedSS)
  {
    fSimSS = TFile::Open(simMergedSS.c_str(), "READ");
    if (fSimSS && !fSimSS->IsZombie())
    {
      simTopSS = fSimSS->GetDirectory(kDirSIM.c_str());
      if (!simTopSS)
      {
        cout << ANSI_BOLD_YEL
             << "[WARN] SS_QA templates: missing topDir '" << kDirSIM
             << "' in merged SIM file: " << simMergedSS
             << ANSI_RESET << "\n";
      }
    }
    else
    {
      cout << ANSI_BOLD_YEL
           << "[WARN] SS_QA templates: cannot open merged SIM file: " << simMergedSS
           << ANSI_RESET << "\n";
      if (fSimSS) { fSimSS->Close(); delete fSimSS; }
      fSimSS = nullptr;
    }
  }
  else
  {
    cout << ANSI_BOLD_YEL
         << "[WARN] SS_QA templates: could not build merged SIM5+10+20 file: " << simMergedSS
         << ANSI_RESET << "\n";
  }
}

if (!skipToCentralityAndPtOverlaysWithSSQA)
{
    TFile* fSimSSEmbedded = nullptr;
    TDirectory* simTopSSEmbedded = nullptr;
    bool attemptedEmbeddedTemplate = false;

    auto CloseTemplateFile = [&](TFile*& f, TDirectory*& d) -> void
    {
        d = nullptr;
        if (f)
        {
            f->Close();
            delete f;
            f = nullptr;
        }
    };

    auto EnsureEmbeddedTemplateTopDir = [&]() -> TDirectory*
    {
        if (attemptedEmbeddedTemplate) return simTopSSEmbedded;
        attemptedEmbeddedTemplate = true;

        if (simTopSSEmbedded) return simTopSSEmbedded;

        const SimSample activeSS = CurrentSimSample();
        if (!IsEmbeddedSimSample(activeSS)) return nullptr;

        string embeddedPath = SimInputPathForSample(activeSS);

        if (activeSS == SimSample::kEmbeddedPhoton10And20Merged)
        {
            if (embeddedPath.empty() || gSystem->AccessPathName(embeddedPath.c_str()))
            {
                const string outMerged =
                    MergedSimEmbeddedPath("photonJet10and20merged_SIM", "RecoilJets_embeddedPhoton10plus20_MERGED.root");

                const bool okMergeEmbedded = BuildMergedSIMFile_PhotonSlices(
                    {InputSimEmbeddedSample("embeddedPhoton10"), InputSimEmbeddedSample("embeddedPhoton20")},
                    {kSigmaPhoton10_pb, kSigmaPhoton20_pb},
                    outMerged,
                    kDirSIM,
                    {"embeddedPhoton10", "embeddedPhoton20"}
                );

                if (okMergeEmbedded && !gSystem->AccessPathName(outMerged.c_str()))
                {
                    embeddedPath = outMerged;
                }
            }
        }

        if (embeddedPath.empty() || gSystem->AccessPathName(embeddedPath.c_str()))
        {
            cout << ANSI_BOLD_YEL
                 << "[WARN] SS_QA embedded templates: active embedded SIM file is unavailable: "
                 << embeddedPath
                 << ANSI_RESET << "\n";
            return nullptr;
        }

        fSimSSEmbedded = TFile::Open(embeddedPath.c_str(), "READ");
        if (!fSimSSEmbedded || fSimSSEmbedded->IsZombie())
        {
            cout << ANSI_BOLD_YEL
                 << "[WARN] SS_QA embedded templates: could not open embedded SIM file: "
                 << embeddedPath
                 << ANSI_RESET << "\n";
            CloseTemplateFile(fSimSSEmbedded, simTopSSEmbedded);
            return nullptr;
        }

        simTopSSEmbedded = fSimSSEmbedded->GetDirectory(kDirSIM.c_str());
        if (!simTopSSEmbedded)
        {
            cout << ANSI_BOLD_YEL
                 << "[WARN] SS_QA embedded templates: missing topDir '" << kDirSIM
                 << "' in embedded SIM file: " << embeddedPath
                 << ANSI_RESET << "\n";
            CloseTemplateFile(fSimSSEmbedded, simTopSSEmbedded);
            return nullptr;
        }

        return simTopSSEmbedded;
    };

    // --- Build vector of all SIM template sources for Pythia overlay tables ---
    struct SSTemplateSource {
        string folderName;
        TDirectory* topDir;
        TFile* ownedFile;   // non-null → we must close it
    };
    vector<SSTemplateSource> ssTemplateSources;

    // 1. Always: merged pp photonJet5+10+20
    if (simTopSS)
        ssTemplateSources.push_back({"overlaysWithPythia_photonJet5and10and20merged", simTopSS, nullptr});

    // 2. Photon+Jet embedded (whichever toggle is active)
    {
        TDirectory* photEmb = EnsureEmbeddedTemplateTopDir();
        if (photEmb)
        {
            string tag = "embeddedPhoton20";
            if (bothPhoton10and20simEmbedded)      tag = "embeddedPhoton10and20merged";
            else if (isPhotonJet10Embedded)         tag = "embeddedPhoton10";
            ssTemplateSources.push_back({"overlaysWithPythia_" + tag, photEmb, nullptr});
        }
    }

    // 3. Inclusive jet embedded (whichever toggle is active)
    TFile* fInclJetSS = nullptr;
    {
        string inclPath;
        string inclTag;
        if (bothInclusiveJet10and20simEmbedded)
        {
            // TODO: merged inclusive jet embedded
        }
        else if (isInclusiveJet20Embedded)
        {
            inclPath = InputInclusiveJetEmbeddedSample("embeddedJet20");
            inclTag  = "embeddedJet20";
        }
        else if (isInclusiveJet10Embedded)
        {
            inclPath = InputInclusiveJetEmbeddedSample("embeddedJet10");
            inclTag  = "embeddedJet10";
        }
        if (!inclPath.empty() && !gSystem->AccessPathName(inclPath.c_str()))
        {
            fInclJetSS = TFile::Open(inclPath.c_str(), "READ");
            if (fInclJetSS && !fInclJetSS->IsZombie())
            {
                TDirectory* d = fInclJetSS->GetDirectory(kDirSIM.c_str());
                if (d)
                    ssTemplateSources.push_back({"overlaysWithPythia_" + inclTag, d, fInclJetSS});
                else
                { fInclJetSS->Close(); delete fInclJetSS; fInclJetSS = nullptr; }
            }
            else
            { if (fInclJetSS) { fInclJetSS->Close(); delete fInclJetSS; } fInclJetSS = nullptr; }
        }
    }

    struct SSOverlayVariantCfg

    {
        string folder;
        vector<std::size_t> indices;
    };

    const vector<SSOverlayVariantCfg> ueOverlayCfgs = {
        {"noSub_baseVariant", {std::size_t(0), std::size_t(1)}},
        {"noSub_variantA",    {std::size_t(0), std::size_t(2)}},
        {"noSub_variantB",    {std::size_t(0), std::size_t(3)}},
        {"variantA_variantB", {std::size_t(2), std::size_t(3)}}
    };

    const vector<SSOverlayVariantCfg> pythiaOverlayCfgs = {
        {"noSub",                 {std::size_t(0)}},
        {"baseVariant",           {std::size_t(1)}},
        {"variantA",              {std::size_t(2)}},
        {"variantB",              {std::size_t(3)}},
        {"noSub_baseVariant",     {std::size_t(0), std::size_t(1)}},
        {"noSub_variantA",        {std::size_t(0), std::size_t(2)}},
        {"variantA_variantB",     {std::size_t(2), std::size_t(3)}},
        {"variantA_variantB_noSub", {std::size_t(0), std::size_t(2), std::size_t(3)}}
    };

    auto HaveAllVariantFiles = [&](const vector<std::size_t>& indices) -> bool
    {
        for (std::size_t idx : indices)
        {
            if (idx >= handles.size()) return false;
            if (!handles[idx].file) return false;
        }
        return true;
    };

    auto DrawUEOverlayTable =
      [&](const SSOverlayVariantCfg& cfg,
          const CentBin& cb,
          const PtBin& b,
          const string& outPng) -> void
    {
        if (!HaveAllVariantFiles(cfg.indices)) return;

        const string cName = TString::Format("c_ssQA_UE_%s_%s_%s_%s",
                                             cfg.folder.c_str(), trigAA.c_str(), cb.folder.c_str(), b.folder.c_str()).Data();

        TCanvas cSS(cName.c_str(), cName.c_str(), 2600, 780);
        cSS.Divide(5, 1, 0.001, 0.001);

        vector<TH1*> keepH;
        vector<TLegend*> keepLeg;
        keepH.reserve(ssVars.size() * cfg.indices.size());
        keepLeg.reserve(ssVars.size());

        bool anyPad = false;

        for (int iv = 0; iv < (int)ssVars.size(); ++iv)
        {
            cSS.cd(iv + 1);
            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.14);
            gPad->SetTopMargin(0.24);
            gPad->SetLogy(false);

            const string hName = "h_ss_" + ssVars[iv].var + "_inclusive" + b.suffix + cb.suffix;

            vector<TH1*> padHists;
            vector<string> padLabels;
            double yMaxPad = 0.0;

            for (std::size_t idx : cfg.indices)
            {
                if (idx >= handles.size()) continue;
                auto& H = handles[idx];
                if (!H.file) continue;

                TDirectory* aaTopSS = H.file->GetDirectory(trigAA.c_str());
                if (!aaTopSS) continue;

                TH1* rawAA = GetTH1FromTopDir(aaTopSS, hName);
                if (!rawAA) continue;

                TH1* hAA = CloneNormalizeStyle(
                    rawAA,
                    TString::Format("ssQA_ueOnly_%s_%s_%s_%s_%s_%zu",
                                    cfg.folder.c_str(), ssVars[iv].var.c_str(), cb.folder.c_str(),
                                    b.folder.c_str(), trigAA.c_str(), idx).Data(),
                    ssColors[idx], ssMarkers[idx]
                );
                if (!hAA) continue;

                hAA->GetXaxis()->UnZoom();
                hAA->SetFillStyle(0);
                hAA->SetLineWidth(2);
                hAA->SetMarkerSize(0.95);

                for (int ib = 1; ib <= hAA->GetNbinsX(); ++ib)
                {
                    yMaxPad = std::max(yMaxPad, (double)(hAA->GetBinContent(ib) + hAA->GetBinError(ib)));
                }

                padHists.push_back(hAA);
                padLabels.push_back(H.label);
                keepH.push_back(hAA);
            }

            if (padHists.empty())
            {
                DrawMissingPad(TString::Format("%s, %d-%d%%, %d-%d GeV",
                                               ssVars[iv].label.c_str(), cb.lo, cb.hi, b.lo, b.hi).Data());
                continue;
            }

            anyPad = true;

            TH1* hFrame = padHists[0];
            hFrame->SetTitle("");
            hFrame->GetXaxis()->SetTitle(ssVars[iv].label.c_str());
            hFrame->GetYaxis()->SetTitle("Normalized");
            hFrame->GetXaxis()->SetTitleSize(0.055);
            hFrame->GetYaxis()->SetTitleSize(0.055);
            hFrame->GetXaxis()->SetLabelSize(0.044);
            hFrame->GetYaxis()->SetLabelSize(0.040);
            hFrame->GetYaxis()->SetTitleOffset(1.55);
            hFrame->SetMinimum(0.0);
            hFrame->SetMaximum((yMaxPad > 0.0) ? (1.12 * yMaxPad) : 1.0);
            hFrame->Draw("E1");
            for (std::size_t ih = 1; ih < padHists.size(); ++ih) padHists[ih]->Draw("E1 same");
            hFrame->Draw("E1 same");

            TLegend* leg = new TLegend(0.08, 0.90, 0.92, 0.985);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.036);
            leg->SetNColumns((int)std::min<std::size_t>(cfg.indices.size(), 2));
            for (std::size_t ih = 0; ih < padHists.size(); ++ih)
            {
                leg->AddEntry(padHists[ih], padLabels[ih].c_str(), "ep");
            }
            leg->Draw();
            keepLeg.push_back(leg);

            TLatex tPad;
            tPad.SetNDC(true);
            tPad.SetTextFont(42);
            tPad.SetTextAlign(22);
            tPad.SetTextSize(0.050);
            tPad.DrawLatex(0.50, 0.845,
                           TString::Format("%s, p_{T}^{#gamma} %d-%d GeV, %d-%d%%",
                                           ssVars[iv].label.c_str(),
                                           b.lo, b.hi, cb.lo, cb.hi).Data());
        }

        if (anyPad)
        {
            SaveCanvas(cSS, outPng);
        }

        for (TLegend* l : keepLeg) delete l;
        for (TH1* h : keepH) delete h;
    };

    auto DrawPythiaOverlaySet =
      [&](const SSOverlayVariantCfg& cfg,
          const CentBin& cb,
          const PtBin& b,
          const string& outRoot,
          TDirectory* templateTopDir) -> void
    {
        if (!HaveAllVariantFiles(cfg.indices)) return;
        if (!ppTop || !templateTopDir) return;

        EnsureDir(outRoot);

        for (const auto& tag : ppg12Tags)
        {
            auto DrawVariantPad =
              [&](bool drawLegend,
                  bool doZoom,
                  int iv,
                  std::vector<TH1*>& keepAlive,
                  std::vector<TLegend*>& keepLegs) -> bool
            {
                const std::string& var = ssVars[iv].var;
                const std::string& vlabel = ssVars[iv].label;
                const bool isW = (var == "weta" || var == "wphi");

                const string hPPName  = string("h_ss_") + var + string("_") + tag + b.suffix;
                const string hSigName = string("h_ss_") + var + string("_") + tag + string("_sig") + b.suffix;
                const string hBkgName = string("h_ss_") + var + string("_") + tag + string("_bkg") + b.suffix;
                const string hAAName  = string("h_ss_") + var + string("_") + tag + b.suffix + cb.suffix;

                TH1* rawPP  = GetTH1FromTopDir(ppTop, hPPName);
                TH1* rawSig = GetTH1FromTopDir(templateTopDir, hSigName);
                TH1* rawBkg = GetTH1FromTopDir(templateTopDir, hBkgName);

                vector<std::pair<std::size_t, TH1*> > rawAAs;
                for (std::size_t idx : cfg.indices)
                {
                    if (idx >= handles.size()) continue;
                    auto& H = handles[idx];
                    if (!H.file) continue;

                    TDirectory* aaTopSS = H.file->GetDirectory(trigAA.c_str());
                    if (!aaTopSS) continue;

                    TH1* rawAA = GetTH1FromTopDir(aaTopSS, hAAName);
                    if (rawAA) rawAAs.push_back(std::make_pair(idx, rawAA));
                }

                if (!rawPP && !rawSig && !rawBkg && rawAAs.empty())
                {
                    DrawMissingPad(TString::Format("%s, %s, %s, %s", var.c_str(), tag.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                    return false;
                }

                TH1* hPP  = nullptr;
                TH1* hSig = nullptr;
                TH1* hBkg = nullptr;
                vector<std::pair<std::size_t, TH1*> > hAAs;

                if (rawPP)
                {
                    hPP = CloneNormalizeStyle(
                        rawPP,
                        TString::Format("ssQA_%s_%s_%s_%s_%s_pp",
                                        cfg.folder.c_str(), tag.c_str(), var.c_str(), b.folder.c_str(), doZoom ? "zoom" : "full").Data(),
                        kRed + 1, 24
                    );

                    if (hPP)
                    {
                        hPP->SetLineWidth(2);
                        hPP->SetLineColor(kRed + 1);
                        hPP->SetMarkerColor(kRed + 1);
                        hPP->SetMarkerStyle(24);
                        hPP->SetMarkerSize(1.00);
                        hPP->SetFillStyle(0);
                    }
                }

                if (rawSig)
                {
                    hSig = CloneNormalizeStyle(
                        rawSig,
                        TString::Format("ssQA_%s_%s_%s_%s_%s_sig",
                                        cfg.folder.c_str(), tag.c_str(), var.c_str(), b.folder.c_str(), doZoom ? "zoom" : "full").Data(),
                        kPink + 7, 24
                    );

                    if (hSig)
                    {
                        hSig->SetLineWidth(2);
                        hSig->SetLineColor(kPink + 7);
                        hSig->SetMarkerColor(kPink + 7);
                        hSig->SetMarkerStyle(1);
                        hSig->SetMarkerSize(0.0);
                        hSig->SetFillStyle(0);
                    }
                }

                if (rawBkg)
                {
                    hBkg = CloneNormalizeStyle(
                        rawBkg,
                        TString::Format("ssQA_%s_%s_%s_%s_%s_bkg",
                                        cfg.folder.c_str(), tag.c_str(), var.c_str(), b.folder.c_str(), doZoom ? "zoom" : "full").Data(),
                        kBlue + 1, 25
                    );

                    if (hBkg)
                    {
                        hBkg->SetLineWidth(2);
                        hBkg->SetLineColor(kBlue + 1);
                        hBkg->SetMarkerColor(kBlue + 1);
                        hBkg->SetMarkerStyle(1);
                        hBkg->SetMarkerSize(0.0);
                        hBkg->SetFillStyle(0);
                    }
                }

                for (const auto& rawPair : rawAAs)
                {
                    const std::size_t idx = rawPair.first;
                    TH1* hAA = CloneNormalizeStyle(
                        rawPair.second,
                        TString::Format("ssQA_%s_%s_%s_%s_%s_aa_%zu",
                                        cfg.folder.c_str(), tag.c_str(), var.c_str(), b.folder.c_str(), doZoom ? "zoom" : "full", idx).Data(),
                        ssColors[idx], ssMarkers[idx]
                    );

                    if (hAA)
                    {
                        hAA->SetLineWidth(2);
                        hAA->SetLineColor(ssColors[idx]);
                        hAA->SetMarkerColor(ssColors[idx]);
                        hAA->SetMarkerStyle(ssMarkers[idx]);
                        hAA->SetMarkerSize(1.00);
                        hAA->SetFillStyle(0);
                        hAAs.push_back(std::make_pair(idx, hAA));
                    }
                }

                TH1* hFrame = nullptr;
                if (doZoom && tag != "tight")
                    hFrame = (hBkg ? hBkg : (!hAAs.empty() ? hAAs.front().second : hPP));
                else
                    hFrame = (hSig ? hSig : (hBkg ? hBkg : (!hAAs.empty() ? hAAs.front().second : hPP)));

                if (!hFrame)
                {
                    DrawMissingPad(TString::Format("%s, %s, %s, %s", var.c_str(), tag.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                    return false;
                }

                hFrame->GetXaxis()->SetTitle(vlabel.c_str());
                hFrame->GetYaxis()->SetTitle("Unit Normalized");
                hFrame->GetYaxis()->SetTitleOffset(1.58);
                hFrame->GetYaxis()->SetTitleSize(0.050);
                hFrame->GetYaxis()->SetLabelSize(0.040);

                double yMax = 0.0;
                if (hPP)
                {
                    for (int ib = 1; ib <= hPP->GetNbinsX(); ++ib)
                        yMax = std::max(yMax, (double)(hPP->GetBinContent(ib) + hPP->GetBinError(ib)));
                }
                if (hSig) yMax = std::max(yMax, (double)hSig->GetMaximum());
                if (hBkg) yMax = std::max(yMax, (double)hBkg->GetMaximum());
                for (const auto& hAAPair : hAAs)
                {
                    TH1* hAA = hAAPair.second;
                    if (!hAA) continue;
                    for (int ib = 1; ib <= hAA->GetNbinsX(); ++ib)
                        yMax = std::max(yMax, (double)(hAA->GetBinContent(ib) + hAA->GetBinError(ib)));
                }

                if (doZoom)
                {
                    double xMinZoom = std::numeric_limits<double>::max();
                    double xMaxZoom = -std::numeric_limits<double>::max();
                    double yMaxZoom = 0.0;

                    auto AccumulateZoomRange = [&](TH1* h) -> void
                    {
                        if (!h) return;
                        for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                        {
                            const double y = h->GetBinContent(ib) + h->GetBinError(ib);
                            if (y <= 0.0) continue;
                            xMinZoom = std::min(xMinZoom, h->GetXaxis()->GetBinLowEdge(ib));
                            xMaxZoom = std::max(xMaxZoom, h->GetXaxis()->GetBinUpEdge(ib));
                            yMaxZoom = std::max(yMaxZoom, y);
                        }
                    };

                    AccumulateZoomRange(hPP);
                    for (const auto& hAAPair : hAAs) AccumulateZoomRange(hAAPair.second);
                    if (tag == "tight") AccumulateZoomRange(hSig);
                    else                  AccumulateZoomRange(hBkg);

                    if (std::isfinite(xMinZoom) && std::isfinite(xMaxZoom) && xMaxZoom > xMinZoom)
                    {
                        const double axisMin = hFrame->GetXaxis()->GetXmin();
                        const double axisMax = hFrame->GetXaxis()->GetXmax();
                        double margin = 0.06 * (xMaxZoom - xMinZoom);
                        if (!(margin > 0.0)) margin = 0.05 * (axisMax - axisMin);

                        const double xLo = std::max(axisMin, xMinZoom - margin);
                        const double xHi = std::min(axisMax, xMaxZoom + margin);
                        hFrame->GetXaxis()->SetRangeUser(xLo, xHi);
                    }

                    if (yMaxZoom > 0.0) yMax = yMaxZoom;
                }

                hFrame->SetMinimum(0.0);
                hFrame->SetMaximum((yMax > 0.0) ? (yMax * (doZoom ? 1.08 : 1.10)) : 1.0);

                if (doZoom)
                {
                    if (tag == "tight")
                    {
                        if (hSig) hSig->Draw("HIST");
                        else if (!hAAs.empty()) hAAs.front().second->Draw("E1");
                        else if (hPP) hPP->Draw("E1");
                    }
                    else
                    {
                        if (hBkg) hBkg->Draw("HIST");
                        else if (!hAAs.empty()) hAAs.front().second->Draw("E1");
                        else if (hPP) hPP->Draw("E1");
                    }
                }
                else
                {
                    if (hSig) hSig->Draw("HIST");
                    else if (hBkg) hBkg->Draw("HIST");
                    else if (!hAAs.empty()) hAAs.front().second->Draw("E1");
                    else if (hPP) hPP->Draw("E1");

                    if (hBkg && hBkg != hSig) hBkg->Draw("HIST same");
                }

                for (const auto& hAAPair : hAAs) hAAPair.second->Draw("E1 same");
                if (hPP) hPP->Draw("E1 same");

                if (drawLegend)
                {
                    TLegend* leg = (isW ? new TLegend(0.41, 0.55, 0.85, 0.80) : new TLegend(0.18, 0.55, 0.62, 0.80));
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(42);
                    leg->SetTextSize(0.036);

                    if (hPP)  leg->AddEntry(hPP,  "pp", "ep");
                    if (hSig) leg->AddEntry(hSig, "Signal MC", "l");
                    if (hBkg) leg->AddEntry(hBkg, "Background MC", "l");
                    for (const auto& hAAPair : hAAs)
                    {
                        leg->AddEntry(
                            hAAPair.second,
                            TString::Format("AuAu (%d-%d%%) %s", cb.lo, cb.hi, handles[hAAPair.first].label.c_str()).Data(),
                            "ep"
                        );
                    }

                    leg->Draw();
                    keepLegs.push_back(leg);
                }

                {
                    TLatex th;
                    th.SetNDC(true);
                    th.SetTextFont(42);
                    if (doZoom)
                    {
                        th.SetTextAlign(33);
                        th.SetTextSize(0.060);
                        th.DrawLatex(0.93, 0.93, "Zoomed");
                    }
                    else
                    {
                        th.SetTextAlign(22);
                        th.SetTextSize(0.050);
                        th.DrawLatex(0.50, 0.91,
                                     TString::Format("%s, %s, p_{T}^{#gamma}: %d-%d GeV",
                                                     vlabel.c_str(), TagLabel(tag).c_str(), b.lo, b.hi).Data());
                    }
                }

                {
                    TLatex tcut;
                    tcut.SetNDC(true);
                    tcut.SetTextFont(42);
                    tcut.SetTextAlign(13);
                    tcut.SetTextSize(doZoom ? 0.034 : 0.040);

                    bool drawCuts = false;
                    bool drawSingleCut = false;
                    double cutLo = 0.0;
                    double cutHi = 0.0;
                    std::string cutText;

                    const bool isPre = (tag == "pre");

                    if (var == "e11e33")
                    {
                        if (isPre)
                        {
                            cutText = "pp presel: #frac{E_{11}}{E_{33}} < 0.98";
                            drawSingleCut = true;
                            cutHi = 0.98;
                        }
                        else
                        {
                            cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                            drawCuts = true;
                            cutLo = 0.4;
                            cutHi = 0.98;
                        }
                    }
                    else if (var == "e32e35")
                    {
                        if (isPre)
                        {
                            cutText = "pp presel: 0.8 < #frac{E_{32}}{E_{35}} < 1.0";
                            drawCuts = true;
                            cutLo = 0.8;
                            cutHi = 1.0;
                        }
                        else
                        {
                            cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                            drawCuts = true;
                            cutLo = 0.92;
                            cutHi = 1.0;
                        }
                    }
                    else if (var == "et1")
                    {
                        if (isPre)
                        {
                            cutText = "pp presel: 0.6 < et1 < 1.0";
                            drawCuts = true;
                            cutLo = 0.6;
                            cutHi = 1.0;
                        }
                        else
                        {
                            cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                            drawCuts = true;
                            cutLo = 0.9;
                            cutHi = 1.0;
                        }
                    }
                    else if (var == "weta")
                    {
                        if (isPre)
                        {
                            cutText = "pp presel: w_{#eta} < 0.6";
                            drawSingleCut = true;
                            cutHi = 0.6;
                        }
                        else
                        {
                            cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                            drawSingleCut = true;
                            const double ptCenter = 0.5 * (b.lo + b.hi);
                            cutHi = 0.15 + 0.006 * ptCenter;
                        }
                    }
                    else if (var == "wphi")
                    {
                        if (!isPre)
                        {
                            cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                            drawSingleCut = true;
                            const double ptCenter = 0.5 * (b.lo + b.hi);
                            cutHi = 0.15 + 0.006 * ptCenter;
                        }
                    }

                    if (!cutText.empty() && !doZoom)
                    {
                        tcut.DrawLatex(0.16, 0.86, cutText.c_str());
                    }

                    if (drawCuts || drawSingleCut)
                    {
                        gPad->Update();
                        const double yMin = gPad->GetUymin();
                        const double yMaxPad = gPad->GetUymax();

                        if (drawCuts)
                        {
                            TLine l1(cutLo, yMin, cutLo, yMaxPad);
                            l1.SetLineColor(kBlack);
                            l1.SetLineWidth(2);
                            l1.SetLineStyle(2);
                            l1.DrawClone("same");

                            TLine l2(cutHi, yMin, cutHi, yMaxPad);
                            l2.SetLineColor(kBlack);
                            l2.SetLineWidth(2);
                            l2.SetLineStyle(2);
                            l2.DrawClone("same");
                        }

                        if (drawSingleCut)
                        {
                            TLine l1(cutHi, yMin, cutHi, yMaxPad);
                            l1.SetLineColor(kBlack);
                            l1.SetLineWidth(2);
                            l1.SetLineStyle(2);
                            l1.DrawClone("same");
                        }
                    }

                    gPad->RedrawAxis();
                }

                if (hPP)  keepAlive.push_back(hPP);
                if (hSig) keepAlive.push_back(hSig);
                if (hBkg) keepAlive.push_back(hBkg);
                for (const auto& hAAPair : hAAs) keepAlive.push_back(hAAPair.second);

                return true;
            };

            TCanvas cPP(
                        TString::Format("c_ssQA_ppDataSigBkg_%s_%s_%s_%s",
                                        cfg.folder.c_str(), tag.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                        "c_ssQA_ppDataSigBkg", 2600, 750
                        );
            cPP.Divide(5, 1, 0.001, 0.001);

            std::vector<TH1*> keepAlive;
            std::vector<TLegend*> keepLegs;
            bool anyPad = false;

            for (int iv = 0; iv < (int)ssVars.size(); ++iv)
            {
                cPP.cd(iv + 1);
                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.18);
                gPad->SetLogy(false);

                if (DrawVariantPad(true, false, iv, keepAlive, keepLegs)) anyPad = true;
            }

            if (anyPad)
            {
                SaveCanvas(cPP, JoinPath(outRoot,
                                         TString::Format("table1x5_PP_SS_%s_DataSigBkg.png", tag.c_str()).Data()));
            }

            for (TLegend* l : keepLegs) delete l;
            keepLegs.clear();
            for (TH1* h : keepAlive) delete h;
            keepAlive.clear();

            TCanvas cPP2(
                         TString::Format("c_ssQA_ppDataSigBkg_zoom_%s_%s_%s_%s",
                                         cfg.folder.c_str(), tag.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                         "c_ssQA_ppDataSigBkg_zoom", 2600, 1200
                         );
            cPP2.Divide(5, 2, 0.001, 0.001);

            std::vector<TH1*> keepAlive2;
            std::vector<TLegend*> keepLegs2;
            bool anyPad2 = false;

            for (int iv = 0; iv < (int)ssVars.size(); ++iv)
            {
                cPP2.cd(iv + 1);
                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.10);
                gPad->SetTopMargin(0.18);
                gPad->SetLogy(false);

                if (DrawVariantPad(true, false, iv, keepAlive2, keepLegs2)) anyPad2 = true;
            }

            for (int iv = 0; iv < (int)ssVars.size(); ++iv)
            {
                cPP2.cd(iv + 6);
                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.06);
                gPad->SetLogy(false);

                if (DrawVariantPad(false, true, iv, keepAlive2, keepLegs2)) anyPad2 = true;
            }

            if (anyPad2)
            {
                SaveCanvas(cPP2, JoinPath(outRoot,
                                          TString::Format("table2x5_PP_SS_%s_DataSigBkg_zoomedPPAuAu.png", tag.c_str()).Data()));
            }

            for (TLegend* l : keepLegs2) delete l;
            keepLegs2.clear();
            for (TH1* h : keepAlive2) delete h;
            keepAlive2.clear();
        }
    };

    for (std::size_t ic = 0; ic < centBins.size(); ++ic)
    {
        const auto& cb = centBins[ic];
        const string centDir = JoinPath(ssQADir, cb.folder);
        EnsureDir(centDir);

        for (int ipt = 0; ipt < kNPtBins; ++ipt)
        {
            const PtBin& b = PtBins()[ipt];
            const string ptDir = JoinPath(centDir, b.folder);
            EnsureDir(ptDir);

            const string ueOverlayDir = JoinPath(ptDir, "UEoverlays");
            EnsureDir(ueOverlayDir);

            for (const auto& cfg : ueOverlayCfgs)
            {
                if (!HaveAllVariantFiles(cfg.indices)) continue;

                const string cfgDir = JoinPath(ueOverlayDir, cfg.folder);
                EnsureDir(cfgDir);
                DrawUEOverlayTable(cfg, cb, b, JoinPath(cfgDir, "table1x5_SS_UEvariantOverlay.png"));
            }

            for (const auto& src : ssTemplateSources)
            {
                if (!src.topDir) continue;
                const string srcDir = JoinPath(ptDir, src.folderName);
                EnsureDir(srcDir);

                for (const auto& cfg : pythiaOverlayCfgs)
                {
                    if (!HaveAllVariantFiles(cfg.indices)) continue;

                    const string cfgDir = JoinPath(srcDir, cfg.folder);
                    EnsureDir(cfgDir);
                    DrawPythiaOverlaySet(cfg, cb, b, cfgDir, src.topDir);
                }
            }
        }
    }

    CloseTemplateFile(fSimSSEmbedded, simTopSSEmbedded);
    if (fInclJetSS) { fInclJetSS->Close(); delete fInclJetSS; fInclJetSS = nullptr; }
}

if (false && !skipToCentralityAndPtOverlaysWithSSQA)
{
    for (std::size_t ic = 0; ic < centBins.size(); ++ic)
    {
        const auto& cb = centBins[ic];
        const string centDir = JoinPath(ssQADir, cb.folder);
        EnsureDir(centDir);
        
        for (int ipt = 0; ipt < kNPtBins; ++ipt)
        {
            const PtBin& b = PtBins()[ipt];
            const string ptDir = JoinPath(centDir, b.folder);
            EnsureDir(ptDir);
            
            const string cName = TString::Format("c_ssQA_1x5_%s_%s_%s",
                                                 trigAA.c_str(), cb.folder.c_str(), b.folder.c_str()).Data();
            
            TCanvas cSS(cName.c_str(), cName.c_str(), 2600, 750);
            cSS.Divide(5, 1, 0.001, 0.001);
            
            vector<TH1*> keepH;
            keepH.reserve(ssVars.size() * 5);
            vector<TLegend*> keepLeg;
            keepLeg.reserve(ssVars.size());
            
            for (int iv = 0; iv < (int)ssVars.size(); ++iv)
            {
                cSS.cd(iv + 1);
                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.16);
                gPad->SetLogy(false);
                
                const string histBase = "h_ss_" + ssVars[iv].var + "_inclusive";
                
                // Collect AuAu variants + PP into per-pad overlays
                vector<TH1*> padHists;
                vector<string> padLabels;
                
                for (std::size_t iu : ssUEIndices)
                {
                    if (iu >= handles.size()) continue;
                    auto& H = handles[iu];
                    if (!H.file) continue;
                    
                    TDirectory* aaTopSS = H.file->GetDirectory(trigAA.c_str());
                    if (!aaTopSS) continue;
                    
                    const string hAAName = histBase + b.suffix + cb.suffix;
                    TH1* rawAA = dynamic_cast<TH1*>(aaTopSS->Get(hAAName.c_str()));
                    if (!rawAA) continue;
                    
                    TH1* hAA = CloneTH1(rawAA,
                                        TString::Format("ssQA_%s_%s_%s_%s_%s",
                                                        ssVars[iv].var.c_str(), H.variant.c_str(),
                                                        cb.folder.c_str(), b.folder.c_str(),
                                                        trigAA.c_str()).Data());
                    if (!hAA) continue;
                    
                    if (ssQA_rebinFactor > 1) hAA->Rebin(ssQA_rebinFactor);
                    EnsureSumw2(hAA);
                    hAA->GetXaxis()->UnZoom();
                    hAA->SetTitle("");
                    
                    const double I = hAA->Integral(0, hAA->GetNbinsX() + 1);
                    if (I > 0.0) hAA->Scale(1.0 / I);
                    
                    hAA->SetLineWidth(2);
                    hAA->SetLineColor(ssColors[iu]);
                    hAA->SetMarkerColor(ssColors[iu]);
                    hAA->SetMarkerStyle(ssMarkers[iu]);
                    hAA->SetMarkerSize(0.95);
                    hAA->SetFillStyle(0);
                    
                    padHists.push_back(hAA);
                    padLabels.push_back(TString::Format("AuAu (%d-%d%%) %s", cb.lo, cb.hi, H.label.c_str()).Data());
                    keepH.push_back(hAA);
                }
                
                // PP overlay
                TH1* hPPss = nullptr;
                if (ppTop)
                {
                    const string hPPName = histBase + b.suffix;
                    TH1* rawPP = dynamic_cast<TH1*>(ppTop->Get(hPPName.c_str()));
                    if (rawPP)
                    {
                        hPPss = CloneTH1(rawPP,
                                         TString::Format("ssQA_pp_%s_%s_%s_%s",
                                                         ssVars[iv].var.c_str(),
                                                         cb.folder.c_str(), b.folder.c_str(),
                                                         trigAA.c_str()).Data());
                        if (hPPss)
                        {
                            if (ssQA_rebinFactor > 1) hPPss->Rebin(ssQA_rebinFactor);
                            EnsureSumw2(hPPss);
                            hPPss->GetXaxis()->UnZoom();
                            hPPss->SetTitle("");
                            
                            const double Ipp = hPPss->Integral(0, hPPss->GetNbinsX() + 1);
                            if (Ipp > 0.0) hPPss->Scale(1.0 / Ipp);
                            
                            hPPss->SetLineWidth(2);
                            hPPss->SetLineColor(kRed + 1);
                            hPPss->SetMarkerColor(kRed + 1);
                            hPPss->SetMarkerStyle(24);
                            hPPss->SetMarkerSize(0.95);
                            hPPss->SetFillStyle(0);
                            
                            padHists.push_back(hPPss);
                            padLabels.push_back("pp");
                            keepH.push_back(hPPss);
                        }
                    }
                }
                
                if (padHists.empty())
                {
                    TLatex tMiss;
                    tMiss.SetNDC(true);
                    tMiss.SetTextFont(42);
                    tMiss.SetTextAlign(22);
                    tMiss.SetTextSize(0.080);
                    tMiss.DrawLatex(0.50, 0.55, "MISSING");
                    continue;
                }
                
                double yMaxPad = 0.0;
                for (auto* h : padHists)
                    yMaxPad = std::max(yMaxPad, h->GetMaximum());
                
                padHists[0]->GetXaxis()->SetTitle(ssVars[iv].label.c_str());
                padHists[0]->GetYaxis()->SetTitle("Normalized");
                padHists[0]->GetYaxis()->SetTitleOffset(1.7);
                padHists[0]->SetMinimum(0.0);
                padHists[0]->SetMaximum(yMaxPad * 1.35);
                
                padHists[0]->Draw("E1");
                for (std::size_t ih = 1; ih < padHists.size(); ++ih)
                    padHists[ih]->Draw("E1 same");
                
                TLegend* leg = new TLegend(0.32, 0.65, 0.85, 0.85);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.045);
                leg->SetNColumns(2);
                for (std::size_t ih = 0; ih < padHists.size(); ++ih)
                    leg->AddEntry(padHists[ih], padLabels[ih].c_str(), "ep");
                leg->Draw();
                keepLeg.push_back(leg);
                
                {
                    TLatex tPad;
                    tPad.SetNDC(true);
                    tPad.SetTextFont(42);
                    tPad.SetTextAlign(22);
                    tPad.SetTextSize(0.058);
                    tPad.DrawLatex(0.50, 0.88,
                                   TString::Format("%s, p_{T}^{#gamma} %d-%d GeV, %d-%d%%",
                                                   ssVars[iv].label.c_str(),
                                                   b.lo, b.hi, cb.lo, cb.hi).Data());
                }
            } // end SS var loop
            
            SaveCanvas(cSS, JoinPath(ptDir, "table1x5_SS_UEvariantOverlay.png"));
            
            for (TLegend* l : keepLeg) delete l;
            for (TH1* h : keepH) delete h;
            
            // PPG12-style tables: pp vs Signal MC vs Background MC + selected AuAu variant(s)
            struct SSOverlayVariantCfg
            {
                string folder;
                vector<std::size_t> indices;
            };
            
            const vector<SSOverlayVariantCfg> ssOverlayVariantCfgs = {
                {"variantA", {std::size_t(2)}},
                {"variantB", {std::size_t(3)}},
                {"variantA_variantB", {std::size_t(2), std::size_t(3)}},
                {"variantA_variantB_noSub", {std::size_t(0), std::size_t(2), std::size_t(3)}}
            };
            
            for (const auto& cfg : ssOverlayVariantCfgs)
            {
                bool haveAnyVariantFile = false;
                for (std::size_t idx : cfg.indices)
                {
                    if (idx < handles.size() && handles[idx].file) { haveAnyVariantFile = true; break; }
                }
                if (!haveAnyVariantFile) continue;
                if (!ppTop || !simTopSS) continue;
                
                const string variantDir = JoinPath(ptDir, cfg.folder);
                EnsureDir(variantDir);
                
                for (const auto& tag : ppg12Tags)
                {
                    auto DrawVariantPad =
                    [&](bool includeTemplates,
                        bool drawLegend,
                        bool doZoom,
                        int iv,
                        std::vector<TH1*>& keepAlive,
                        std::vector<TLegend*>& keepLegPP) -> bool
                    {
                        const std::string& var = ssVars[iv].var;
                        const std::string& vlabel = ssVars[iv].label;
                        
                        const bool isW = (var == "weta" || var == "wphi");
                        
                        const string hPPName  = string("h_ss_") + var + string("_") + tag + b.suffix;
                        const string hSigName = string("h_ss_") + var + string("_") + tag + string("_sig") + b.suffix;
                        const string hBkgName = string("h_ss_") + var + string("_") + tag + string("_bkg") + b.suffix;
                        const string hAAName  = string("h_ss_") + var + string("_") + tag + b.suffix + cb.suffix;
                        
                        TH1* rawPP  = GetTH1FromTopDir(ppTop, hPPName);
                        TH1* rawSig = includeTemplates ? GetTH1FromTopDir(simTopSS, hSigName) : nullptr;
                        TH1* rawBkg = includeTemplates ? GetTH1FromTopDir(simTopSS, hBkgName) : nullptr;
                        
                        vector<std::pair<std::size_t, TH1*> > rawAAs;
                        for (std::size_t idx : cfg.indices)
                        {
                            if (idx >= handles.size()) continue;
                            auto& H = handles[idx];
                            if (!H.file) continue;
                            
                            TDirectory* aaTopSS = H.file->GetDirectory(trigAA.c_str());
                            if (!aaTopSS) continue;
                            
                            TH1* rawAA = GetTH1FromTopDir(aaTopSS, hAAName);
                            if (rawAA) rawAAs.push_back(std::make_pair(idx, rawAA));
                        }
                        
                        if (!rawPP && !rawSig && !rawBkg && rawAAs.empty())
                        {
                            DrawMissingPad(TString::Format("%s, %s, %s, %s", var.c_str(), tag.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                            return false;
                        }
                        
                        TH1* hPP  = nullptr;
                        TH1* hSig = nullptr;
                        TH1* hBkg = nullptr;
                        vector<std::pair<std::size_t, TH1*> > hAAs;
                        
                        if (rawPP)
                        {
                            hPP = CloneNormalizeStyle(rawPP,
                                                      TString::Format("ssQA_%s_%s_%s_%s_%s_pp", cfg.folder.c_str(), tag.c_str(), var.c_str(), b.folder.c_str(), doZoom ? "zoom" : "full").Data(),
                                                      kRed + 1, 24);
                            
                            if (hPP)
                            {
                                hPP->SetLineWidth(2);
                                hPP->SetLineColor(kRed + 1);
                                hPP->SetMarkerColor(kRed + 1);
                                hPP->SetMarkerStyle(24);
                                hPP->SetMarkerSize(1.00);
                                hPP->SetFillStyle(0);
                            }
                        }
                        
                        if (rawSig)
                        {
                            hSig = CloneNormalizeStyle(rawSig,
                                                       TString::Format("ssQA_%s_%s_%s_%s_%s_sig", cfg.folder.c_str(), tag.c_str(), var.c_str(), b.folder.c_str(), doZoom ? "zoom" : "full").Data(),
                                                       kPink + 7, 24);
                            
                            if (hSig)
                            {
                                hSig->SetLineWidth(2);
                                hSig->SetLineColor(kPink + 7);
                                hSig->SetMarkerColor(kPink + 7);
                                hSig->SetMarkerStyle(1);
                                hSig->SetMarkerSize(0.0);
                                hSig->SetFillStyle(0);
                            }
                        }
                        
                        if (rawBkg)
                        {
                            hBkg = CloneNormalizeStyle(rawBkg,
                                                       TString::Format("ssQA_%s_%s_%s_%s_%s_bkg", cfg.folder.c_str(), tag.c_str(), var.c_str(), b.folder.c_str(), doZoom ? "zoom" : "full").Data(),
                                                       kBlue + 1, 25);
                            
                            if (hBkg)
                            {
                                hBkg->SetLineWidth(2);
                                hBkg->SetLineColor(kBlue + 1);
                                hBkg->SetMarkerColor(kBlue + 1);
                                hBkg->SetMarkerStyle(1);
                                hBkg->SetMarkerSize(0.0);
                                hBkg->SetFillStyle(0);
                            }
                        }
                        
                        for (const auto& rawPair : rawAAs)
                        {
                            const std::size_t idx = rawPair.first;
                            TH1* hAA = CloneNormalizeStyle(rawPair.second,
                                                           TString::Format("ssQA_%s_%s_%s_%s_%s_aa_%zu", cfg.folder.c_str(), tag.c_str(), var.c_str(), b.folder.c_str(), doZoom ? "zoom" : "full", idx).Data(),
                                                           ssColors[idx], ssMarkers[idx]);
                            
                            if (hAA)
                            {
                                hAA->SetLineWidth(2);
                                hAA->SetLineColor(ssColors[idx]);
                                hAA->SetMarkerColor(ssColors[idx]);
                                hAA->SetMarkerStyle(ssMarkers[idx]);
                                hAA->SetMarkerSize(1.00);
                                hAA->SetFillStyle(0);
                                hAAs.push_back(std::make_pair(idx, hAA));
                            }
                        }
                        
                        TH1* hFrame = nullptr;
                        if (includeTemplates)
                        {
                            if (doZoom && tag != "tight")
                                hFrame = (hBkg ? hBkg : (!hAAs.empty() ? hAAs.front().second : hPP));
                            else
                                hFrame = (hSig ? hSig : (hBkg ? hBkg : (!hAAs.empty() ? hAAs.front().second : hPP)));
                        }
                        else
                            hFrame = (!hAAs.empty() ? hAAs.front().second : hPP);
                        
                        if (!hFrame)
                        {
                            DrawMissingPad(TString::Format("%s, %s, %s, %s", var.c_str(), tag.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                            return false;
                        }
                        
                        hFrame->GetXaxis()->SetTitle(vlabel.c_str());
                        hFrame->GetYaxis()->SetTitle("Unit Normalized");
                        hFrame->GetYaxis()->SetTitleOffset(1.58);
                        hFrame->GetYaxis()->SetTitleSize(0.050);
                        hFrame->GetYaxis()->SetLabelSize(0.040);
                        
                        double yMax = 0.0;
                        if (includeTemplates)
                        {
                            if (hPP)
                            {
                                for (int ib = 1; ib <= hPP->GetNbinsX(); ++ib)
                                    yMax = std::max(yMax, (double)(hPP->GetBinContent(ib) + hPP->GetBinError(ib)));
                            }
                            if (hSig) yMax = std::max(yMax, (double)hSig->GetMaximum());
                            if (hBkg) yMax = std::max(yMax, (double)hBkg->GetMaximum());
                            for (const auto& hAAPair : hAAs)
                            {
                                TH1* hAA = hAAPair.second;
                                if (!hAA) continue;
                                for (int ib = 1; ib <= hAA->GetNbinsX(); ++ib)
                                    yMax = std::max(yMax, (double)(hAA->GetBinContent(ib) + hAA->GetBinError(ib)));
                            }
                        }
                        else
                        {
                            if (hPP)
                            {
                                for (int ib = 1; ib <= hPP->GetNbinsX(); ++ib)
                                    yMax = std::max(yMax, (double)(hPP->GetBinContent(ib) + hPP->GetBinError(ib)));
                            }
                            for (const auto& hAAPair : hAAs)
                            {
                                TH1* hAA = hAAPair.second;
                                if (!hAA) continue;
                                for (int ib = 1; ib <= hAA->GetNbinsX(); ++ib)
                                    yMax = std::max(yMax, (double)(hAA->GetBinContent(ib) + hAA->GetBinError(ib)));
                            }
                        }
                        
                        if (doZoom)
                        {
                            double xMinZoom = std::numeric_limits<double>::max();
                            double xMaxZoom = -std::numeric_limits<double>::max();
                            double yMaxZoom = 0.0;
                            
                            auto AccumulateZoomRange = [&](TH1* h)
                            {
                                if (!h) return;
                                for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                                {
                                    const double y = h->GetBinContent(ib) + h->GetBinError(ib);
                                    if (y <= 0.0) continue;
                                    xMinZoom = std::min(xMinZoom, h->GetXaxis()->GetBinLowEdge(ib));
                                    xMaxZoom = std::max(xMaxZoom, h->GetXaxis()->GetBinUpEdge(ib));
                                    yMaxZoom = std::max(yMaxZoom, y);
                                }
                            };
                            
                            AccumulateZoomRange(hPP);
                            for (const auto& hAAPair : hAAs) AccumulateZoomRange(hAAPair.second);
                            if (includeTemplates)
                            {
                                if (tag == "tight") AccumulateZoomRange(hSig);
                                else                AccumulateZoomRange(hBkg);
                            }
                            
                            if (std::isfinite(xMinZoom) && std::isfinite(xMaxZoom) && xMaxZoom > xMinZoom)
                            {
                                const double axisMin = hFrame->GetXaxis()->GetXmin();
                                const double axisMax = hFrame->GetXaxis()->GetXmax();
                                double margin = 0.06 * (xMaxZoom - xMinZoom);
                                if (!(margin > 0.0)) margin = 0.05 * (axisMax - axisMin);
                                
                                const double xLo = std::max(axisMin, xMinZoom - margin);
                                const double xHi = std::min(axisMax, xMaxZoom + margin);
                                hFrame->GetXaxis()->SetRangeUser(xLo, xHi);
                            }
                            
                            if (yMaxZoom > 0.0) yMax = yMaxZoom;
                        }
                        
                        hFrame->SetMinimum(0.0);
                        hFrame->SetMaximum((yMax > 0.0) ? (yMax * (doZoom ? 1.08 : 1.10)) : 1.0);
                        
                        if (includeTemplates)
                        {
                            if (doZoom)
                            {
                                // Zoomed row: only the relevant MC template
                                if (tag == "tight")
                                {
                                    if (hSig)             hSig->Draw("HIST");
                                    else if (!hAAs.empty()) hAAs.front().second->Draw("E1");
                                    else if (hPP)         hPP->Draw("E1");
                                }
                                else
                                {
                                    if (hBkg)             hBkg->Draw("HIST");
                                    else if (!hAAs.empty()) hAAs.front().second->Draw("E1");
                                    else if (hPP)         hPP->Draw("E1");
                                }
                            }
                            else
                            {
                                // Full row: both MC templates
                                if (hSig)
                                {
                                    hSig->Draw("HIST");
                                }
                                else if (hBkg)
                                {
                                    hBkg->Draw("HIST");
                                }
                                else if (!hAAs.empty())
                                {
                                    hAAs.front().second->Draw("E1");
                                }
                                else
                                {
                                    hPP->Draw("E1");
                                }
                                
                                if (hBkg && hBkg != hSig) hBkg->Draw("HIST same");
                            }
                        }
                        else
                        {
                            if (!hAAs.empty())
                            {
                                hAAs.front().second->Draw("E1");
                            }
                            else if (hPP)
                            {
                                hPP->Draw("E1");
                            }
                        }
                        
                        for (const auto& hAAPair : hAAs) hAAPair.second->Draw("E1 same");
                        if (hPP) hPP->Draw("E1 same");
                        
                        if (drawLegend)
                        {
                            TLegend* leg = (isW ? new TLegend(0.41, 0.55, 0.85, 0.8) : new TLegend(0.18, 0.55, 0.62, 0.8));
                            leg->SetBorderSize(0);
                            leg->SetFillStyle(0);
                            leg->SetTextFont(42);
                            leg->SetTextSize(0.036);
                            
                            if (hPP)  leg->AddEntry(hPP,  "pp", "ep");
                            if (hSig) leg->AddEntry(hSig, "Signal MC", "l");
                            if (hBkg) leg->AddEntry(hBkg, "Background MC", "l");
                            for (const auto& hAAPair : hAAs)
                            {
                                leg->AddEntry(
                                              hAAPair.second,
                                              TString::Format("AuAu (%d-%d%%) %s", cb.lo, cb.hi, handles[hAAPair.first].label.c_str()).Data(),
                                              "ep"
                                              );
                            }
                            
                            leg->Draw();
                            keepLegPP.push_back(leg);
                        }
                        
                        {
                            TLatex th;
                            th.SetNDC(true);
                            th.SetTextFont(42);
                            if (doZoom)
                            {
                                th.SetTextAlign(33);
                                th.SetTextSize(0.060);
                                th.DrawLatex(0.93, 0.93, "Zoomed");
                            }
                            else
                            {
                                th.SetTextAlign(22);
                                th.SetTextSize(0.050);
                                th.DrawLatex(0.50, 0.91,
                                             TString::Format("%s, %s, p_{T}^{#gamma}: %d-%d GeV",
                                                             vlabel.c_str(), TagLabel(tag).c_str(), b.lo, b.hi).Data());
                            }
                        }
                        
                        {
                            TLatex tcut;
                            tcut.SetNDC(true);
                            tcut.SetTextFont(42);
                            tcut.SetTextAlign(13);
                            tcut.SetTextSize(doZoom ? 0.034 : 0.040);
                            
                            bool drawCuts = false;
                            bool drawSingleCut = false;
                            double cutLo = 0.0;
                            double cutHi = 0.0;
                            std::string cutText;
                            
                            const bool isPre = (tag == "pre");
                            
                            if (var == "e11e33")
                            {
                                if (isPre)
                                {
                                    cutText = "pp presel: #frac{E_{11}}{E_{33}} < 0.98";
                                    drawSingleCut = true;
                                    cutHi = 0.98;
                                }
                                else
                                {
                                    cutText = "Tight #gamma-ID: 0.4 < #frac{E_{11}}{E_{33}} < 0.98";
                                    drawCuts = true;
                                    cutLo = 0.4;
                                    cutHi = 0.98;
                                }
                            }
                            else if (var == "e32e35")
                            {
                                if (isPre)
                                {
                                    cutText = "pp presel: 0.8 < #frac{E_{32}}{E_{35}} < 1.0";
                                    drawCuts = true;
                                    cutLo = 0.8;
                                    cutHi = 1.0;
                                }
                                else
                                {
                                    cutText = "#gamma-ID: 0.92 < #frac{E_{32}}{E_{35}} < 1.0";
                                    drawCuts = true;
                                    cutLo = 0.92;
                                    cutHi = 1.0;
                                }
                            }
                            else if (var == "et1")
                            {
                                if (isPre)
                                {
                                    cutText = "pp presel: 0.6 < et1 < 1.0";
                                    drawCuts = true;
                                    cutLo = 0.6;
                                    cutHi = 1.0;
                                }
                                else
                                {
                                    cutText = "#gamma-ID: 0.9 < et1 < 1.0";
                                    drawCuts = true;
                                    cutLo = 0.9;
                                    cutHi = 1.0;
                                }
                            }
                            else if (var == "weta")
                            {
                                if (isPre)
                                {
                                    cutText = "pp presel: w_{#eta} < 0.6";
                                    drawSingleCut = true;
                                    cutHi = 0.6;
                                }
                                else
                                {
                                    cutText = "#gamma-ID: 0 < w_{#eta}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                                    drawSingleCut = true;
                                    const double ptCenter = 0.5 * (b.lo + b.hi);
                                    cutHi = 0.15 + 0.006 * ptCenter;
                                }
                            }
                            else if (var == "wphi")
                            {
                                if (!isPre)
                                {
                                    cutText = "#gamma-ID: 0 < w_{#phi}^{cogX} < 0.15 + 0.006 E_{T}^{#gamma}";
                                    drawSingleCut = true;
                                    const double ptCenter = 0.5 * (b.lo + b.hi);
                                    cutHi = 0.15 + 0.006 * ptCenter;
                                }
                                // isPre: no preselection cut on wphi — no line, no text
                            }
                            
                            if (!cutText.empty() && !doZoom)
                            {
                                tcut.DrawLatex(0.16, 0.86, cutText.c_str());
                            }
                            
                            if (drawCuts || drawSingleCut)
                            {
                                gPad->Update();
                                const double yMin = gPad->GetUymin();
                                const double yMaxPad = gPad->GetUymax();
                                
                                if (drawCuts)
                                {
                                    TLine l1(cutLo, yMin, cutLo, yMaxPad);
                                    l1.SetLineColor(kBlack);
                                    l1.SetLineWidth(2);
                                    l1.SetLineStyle(2);
                                    l1.DrawClone("same");
                                    
                                    TLine l2(cutHi, yMin, cutHi, yMaxPad);
                                    l2.SetLineColor(kBlack);
                                    l2.SetLineWidth(2);
                                    l2.SetLineStyle(2);
                                    l2.DrawClone("same");
                                }
                                
                                if (drawSingleCut)
                                {
                                    TLine l1(cutHi, yMin, cutHi, yMaxPad);
                                    l1.SetLineColor(kBlack);
                                    l1.SetLineWidth(2);
                                    l1.SetLineStyle(2);
                                    l1.DrawClone("same");
                                }
                            }
                            
                            gPad->RedrawAxis();
                        }
                        
                        if (hPP)  keepAlive.push_back(hPP);
                        if (hSig) keepAlive.push_back(hSig);
                        if (hBkg) keepAlive.push_back(hBkg);
                        for (const auto& hAAPair : hAAs) keepAlive.push_back(hAAPair.second);
                        
                        return true;
                    };
                    
                    TCanvas cPP(
                                TString::Format("c_ssQA_ppDataSigBkg_%s_%s_%s_%s",
                                                cfg.folder.c_str(), tag.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                "c_ssQA_ppDataSigBkg", 2600, 750
                                );
                    cPP.Divide(5, 1, 0.001, 0.001);
                    
                    std::vector<TH1*> keepAlive;
                    keepAlive.reserve(ssVars.size() * 6);
                    
                    std::vector<TLegend*> keepLegPP;
                    keepLegPP.reserve(ssVars.size());
                    
                    bool anyPad = false;
                    
                    for (int iv = 0; iv < (int)ssVars.size(); ++iv)
                    {
                        cPP.cd(iv + 1);
                        gPad->SetLeftMargin(0.14);
                        gPad->SetRightMargin(0.05);
                        gPad->SetBottomMargin(0.14);
                        gPad->SetTopMargin(0.18);
                        gPad->SetLogy(false);
                        
                        if (DrawVariantPad(true, true, false, iv, keepAlive, keepLegPP)) anyPad = true;
                    }
                    
                    if (anyPad)
                    {
                        SaveCanvas(cPP, JoinPath(variantDir,
                                                 TString::Format("table1x5_PP_SS_%s_DataSigBkg.png", tag.c_str()).Data()));
                    }
                    
                    for (TLegend* l : keepLegPP) delete l;
                    keepLegPP.clear();
                    
                    for (TH1* h : keepAlive) delete h;
                    keepAlive.clear();
                    
                    TCanvas cPP2(
                                 TString::Format("c_ssQA_ppDataSigBkg_zoom_%s_%s_%s_%s",
                                                 cfg.folder.c_str(), tag.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                 "c_ssQA_ppDataSigBkg_zoom", 2600, 1200
                                 );
                    cPP2.Divide(5, 2, 0.001, 0.001);
                    std::vector<TH1*> keepAlive2;
                    keepAlive2.reserve(ssVars.size() * 10);
                    
                    std::vector<TLegend*> keepLegPP2;
                    keepLegPP2.reserve(ssVars.size());
                    
                    bool anyPad2 = false;
                    
                    for (int iv = 0; iv < (int)ssVars.size(); ++iv)
                    {
                        cPP2.cd(iv + 1);
                        gPad->SetLeftMargin(0.14);
                        gPad->SetRightMargin(0.05);
                        gPad->SetBottomMargin(0.10);
                        gPad->SetTopMargin(0.18);
                        gPad->SetLogy(false);
                        
                        if (DrawVariantPad(true, true, false, iv, keepAlive2, keepLegPP2)) anyPad2 = true;
                    }
                    
                    for (int iv = 0; iv < (int)ssVars.size(); ++iv)
                    {
                        cPP2.cd(iv + 6);
                        gPad->SetLeftMargin(0.14);
                        gPad->SetRightMargin(0.05);
                        gPad->SetBottomMargin(0.14);
                        gPad->SetTopMargin(0.06);
                        
                        gPad->SetLogy(false);
                        
                        if (DrawVariantPad(true, false, true, iv, keepAlive2, keepLegPP2)) anyPad2 = true;
                    }
                    
                    if (anyPad2)
                    {
                        SaveCanvas(cPP2, JoinPath(variantDir,
                                                  TString::Format("table2x5_PP_SS_%s_DataSigBkg_zoomedPPAuAu.png", tag.c_str()).Data()));
                    }
                    
                    for (TLegend* l : keepLegPP2) delete l;
                    keepLegPP2.clear();
                    
                    for (TH1* h : keepAlive2) delete h;
                    keepAlive2.clear();
                }
            }
        } // end pT loop
    } // end cent loop
}

{
    const string perPtBinOverlayBase = JoinPath(ssQADir, "perPtBinOverlays");
    const string perCentralityOverlayBase = JoinPath(ssQADir, "perCentralityOverlays");
    EnsureDir(perPtBinOverlayBase);
    EnsureDir(perCentralityOverlayBase);
    
    const vector<std::size_t> ssTableVariantIdx = {std::size_t(0), std::size_t(2), std::size_t(3)};
    const int overlayColors[] = {
        kBlack, kOrange + 7, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
        kCyan + 1, kViolet + 1, kAzure + 1, kSpring + 5, kPink + 1, kTeal + 2
    };
    const vector<string> centOverlayTags = {"inclusive", "pre", "tight", "nonTight"};
    
    // --- Open embedded sim files for per-centrality overlay Signal/Bkg MC ---
    struct EmbeddedSSOverlaySource {
        string folderName;
        string label;
        TFile* file = nullptr;
        TDirectory* topDir = nullptr;
    };
    
    auto OpenEmbeddedSSFile = [&](const string& path) -> std::pair<TFile*, TDirectory*>
    {
        if (path.empty() || gSystem->AccessPathName(path.c_str())) return {nullptr, nullptr};
        TFile* f = TFile::Open(path.c_str(), "READ");
        if (!f || f->IsZombie()) { if (f) { f->Close(); delete f; } return {nullptr, nullptr}; }
        TDirectory* d = f->GetDirectory(kDirSIM.c_str());
        if (!d) { f->Close(); delete f; return {nullptr, nullptr}; }
        return {f, d};
    };
    
    vector<EmbeddedSSOverlaySource> embeddedSSSources;
    
    {
        // PhotonJet embedded
        string photonEmbPath;
        string photonEmbLabel;
        if (bothPhoton10and20simEmbedded)
        {
            photonEmbPath = MergedSimEmbeddedPath("photonJet10and20merged_SIM",
                                                   "RecoilJets_embeddedPhoton10plus20_MERGED.root");
            photonEmbLabel = "Photon+Jet Embedded (10+20)";
        }
        else if (isPhotonJet20Embedded)
        {
            photonEmbPath = InputSimEmbeddedSample("embeddedPhoton20");
            photonEmbLabel = "Photon+Jet 20 Embedded";
        }
        else if (isPhotonJet10Embedded)
        {
            photonEmbPath = InputSimEmbeddedSample("embeddedPhoton10");
            photonEmbLabel = "Photon+Jet 10 Embedded";
        }
        if (!photonEmbPath.empty())
        {
            auto [f, d] = OpenEmbeddedSSFile(photonEmbPath);
            if (d) embeddedSSSources.push_back({"PhotonJetEmbedded", photonEmbLabel, f, d});
            else if (f) { f->Close(); delete f; }
        }
        
        // InclusiveJet embedded
        string inclEmbPath;
        string inclEmbLabel;
        if (bothInclusiveJet10and20simEmbedded)
        {
            // TODO: merged inclusive jet embedded if you ever build it
            inclEmbLabel = "Inclusive Jet Embedded (10+20)";
        }
        else if (isInclusiveJet20Embedded)
        {
            inclEmbPath = InputInclusiveJetEmbeddedSample("embeddedJet20");
            inclEmbLabel = "Inclusive Jet 20 Embedded";
        }
        else if (isInclusiveJet10Embedded)
        {
            inclEmbPath = InputInclusiveJetEmbeddedSample("embeddedJet10");
            inclEmbLabel = "Inclusive Jet 10 Embedded";
        }
        if (!inclEmbPath.empty())
        {
            auto [f, d] = OpenEmbeddedSSFile(inclEmbPath);
            if (d) embeddedSSSources.push_back({"inclusiveJetEmbedded", inclEmbLabel, f, d});
            else if (f) { f->Close(); delete f; }
        }
    }
    
    auto GetSSOverlaySelectionText =
    [&](const string& tag,
        bool overlayPtBins,
        double ptCenter) -> string
    {
        if (tag == "pre")
        {
            return "#gamma-presel: #frac{E_{11}}{E_{33}} < 0.98, 0.6 < et1 < 1.0, 0.8 < #frac{E_{32}}{E_{35}} < 1.0, w_{#eta} < 0.6";
        }
        
        if (tag == "tight")
        {
            if (overlayPtBins)
            {
                return "#gamma-tight: 0.4 < #frac{E_{11}}{E_{33}} < 0.98, 0.9 < et1 < 1.0, 0.92 < #frac{E_{32}}{E_{35}} < 1.0, 0 < w_{#eta/#phi} < 0.15 + 0.006 E_{T}^{#gamma}";
            }
            
            const double tightWCut = 0.15 + 0.006 * ptCenter;
            return TString::Format("#gamma-tight: 0.4 < #frac{E_{11}}{E_{33}} < 0.98, 0.9 < et1 < 1.0, 0.92 < #frac{E_{32}}{E_{35}} < 1.0, 0 < w_{#eta/#phi} < 0.15 + 0.006 E_{T}^{#gamma} = %.3f",
                                   tightWCut).Data();
        }
        
        if (tag == "nonTight")
        {
            return "#gamma-nonTight: fail #geq 2 of 5 tight cuts";
        }
        
        return "";
    };
    
    auto DrawSSOverlaySelectionBlock =
    [&](const string& tag,
        bool overlayPtBins,
        double ptCenter) -> bool
    {
        if (overlayPtBins) return false;
        
        TLatex tSel;
        tSel.SetNDC(true);
        tSel.SetTextFont(42);
        tSel.SetTextAlign(13);
        
        const double xLabel = 0.58;
        const double xColL  = 0.58;
        const double xColR  = 0.77;
        const double yLabel = 0.31;
        const double yRow1  = 0.25;
        const double yRow2  = 0.19;
        
        if (tag == "pre")
        {
            tSel.SetTextSize(0.040);
            tSel.DrawLatex(xLabel, yLabel, "#gamma-presel:");
            
            tSel.SetTextSize(0.032);
            tSel.DrawLatex(xColL, yRow1, "#frac{E_{11}}{E_{33}} < 0.98");
            tSel.DrawLatex(xColR, yRow1, "0.6 < et1 < 1.0");
            tSel.DrawLatex(xColL, yRow2, "0.8 < #frac{E_{32}}{E_{35}} < 1.0");
            tSel.DrawLatex(xColR, yRow2, "w_{#eta} < 0.6");
            return true;
        }
        
        if (tag == "tight")
        {
            tSel.SetTextSize(0.040);
            tSel.DrawLatex(xLabel, yLabel, "#gamma-tight:");
            
            tSel.SetTextSize(0.032);
            tSel.DrawLatex(xColL, yRow1, "0.4 < #frac{E_{11}}{E_{33}} < 0.98");
            tSel.DrawLatex(xColR, yRow1, "0.9 < et1 < 1.0");
            tSel.DrawLatex(xColL, yRow2, "0.92 < #frac{E_{32}}{E_{35}} < 1.0");
            tSel.DrawLatex(
                xColR, yRow2,
                TString::Format("0 < w_{#eta/#phi} < 0.15 + 0.006 E_{T}^{#gamma} = %.3f",
                                0.15 + 0.006 * ptCenter).Data()
            );
            return true;
        }
        
        if (tag == "nonTight")
        {
            const double xNT = 0.38;
            tSel.SetTextSize(0.040);
            tSel.DrawLatex(xNT, yLabel, "#gamma-nonTight:");
            
            tSel.SetTextSize(0.032);
            tSel.DrawLatex(xNT, yRow1, "fail #geq 2 of 5 tight cuts");
            return true;
        }
        
        return false;
    };
    
    auto DrawSSOverlayCutsAndText =
    [&](const string& var,
        const string& tag,
        bool overlayPtBins,
        double ptCenter,
        double textX,
        double textY,
        double textSize) -> void
    {
        const bool drewBottomRightBlock = DrawSSOverlaySelectionBlock(tag, overlayPtBins, ptCenter);
        
        const string text = GetSSOverlaySelectionText(tag, overlayPtBins, ptCenter);
        if (!text.empty() && !drewBottomRightBlock)
        {
            double drawTextX = textX;
            double drawTextY = textY;
            double drawTextSize = textSize;
            
            TLatex tSel;
            tSel.SetNDC(true);
            tSel.SetTextFont(42);
            tSel.SetTextAlign(13);
            tSel.SetTextSize(drawTextSize);
            tSel.DrawLatex(drawTextX, drawTextY, text.c_str());
        }
        
        string cutTag = tag;
        if (cutTag == "nonTight") cutTag = "tight";
        if (cutTag != "pre" && cutTag != "tight") return;
        
        bool drawLo = false;
        bool drawHi = false;
        double cutLo = 0.0;
        double cutHi = 0.0;
        
        if (cutTag == "pre")
        {
            if (var == "e11e33")
            {
                drawHi = true;
                cutHi = 0.98;
            }
            else if (var == "et1")
            {
                drawLo = true;
                drawHi = true;
                cutLo = 0.6;
                cutHi = 1.0;
            }
            else if (var == "e32e35")
            {
                drawLo = true;
                drawHi = true;
                cutLo = 0.8;
                cutHi = 1.0;
            }
            else if (var == "weta")
            {
                drawHi = true;
                cutHi = 0.6;
            }
        }
        else if (cutTag == "tight")
        {
            if (var == "e11e33")
            {
                drawLo = true;
                drawHi = true;
                cutLo = 0.4;
                cutHi = 0.98;
            }
            else if (var == "et1")
            {
                drawLo = true;
                drawHi = true;
                cutLo = 0.9;
                cutHi = 1.0;
            }
            else if (var == "e32e35")
            {
                drawLo = true;
                drawHi = true;
                cutLo = 0.92;
                cutHi = 1.0;
            }
            else if (var == "weta" || var == "wphi")
            {
                if (!(overlayPtBins && (var == "weta" || var == "wphi")))
                {
                    drawLo = true;
                    drawHi = true;
                    cutLo = 0.0;
                    cutHi = 0.15 + 0.006 * ptCenter;
                }
            }
        }
        
        if (!drawLo && !drawHi) return;
        
        gPad->Update();
        const double yMin = gPad->GetUymin();
        const double yMax = gPad->GetUymax();
        
        if (drawLo)
        {
            TLine lLo(cutLo, yMin, cutLo, yMax);
            lLo.SetLineColor(kBlack);
            lLo.SetLineWidth(2);
            lLo.SetLineStyle(2);
            lLo.DrawClone("same");
        }
        
        if (drawHi)
        {
            TLine lHi(cutHi, yMin, cutHi, yMax);
            lHi.SetLineColor(kBlack);
            lHi.SetLineWidth(2);
            lHi.SetLineStyle(2);
            lHi.DrawClone("same");
        }
        
        gPad->RedrawAxis();
    };
    
    auto DrawSSOverlayTable3x5 =
    [&](const string& outDir,
        const string& outName,
        const string& tag,
        bool overlayPtBins,
        const CentBin* fixedCent,
        const PtBin* fixedPt) -> void
    {
        TCanvas cTbl(
                     TString::Format("c_ssQA_3x5_%s_%s_%s_%s",
                                     overlayPtBins ? "perPt" : "perCent",
                                     tag.c_str(),
                                     overlayPtBins ? fixedCent->folder.c_str() : fixedPt->folder.c_str(),
                                     trigAA.c_str()).Data(),
                     "c_ssQA_3x5", 2600, 1500
                     );
        cTbl.Divide(5, 3, 0.001, 0.001);
        
        vector<TH1*> keepAlive;
        keepAlive.reserve(ssVars.size() * 3 * 12);
        
        vector<TLegend*> keepLegs;
        keepLegs.reserve(3);
        
        bool anyPad = false;
        const int nOverlayColors = (int)(sizeof(overlayColors) / sizeof(overlayColors[0]));
        
        for (int irow = 0; irow < 3; ++irow)
        {
            const std::size_t vidx = ssTableVariantIdx[(std::size_t)irow];
            if (vidx >= handles.size()) continue;
            auto& H = handles[vidx];
            if (!H.file) continue;
            
            TDirectory* aaTopSS = H.file->GetDirectory(trigAA.c_str());
            if (!aaTopSS) continue;
            
            for (int ivar = 0; ivar < (int)ssVars.size(); ++ivar)
            {
                cTbl.cd(irow * 5 + ivar + 1);
                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.05);
                gPad->SetBottomMargin(0.14);
                gPad->SetTopMargin(0.18);
                gPad->SetLogy(false);
                
                const std::string& var = ssVars[ivar].var;
                const std::string& vlabel = ssVars[ivar].label;
                
                vector<TH1*> hOverlays;
                vector<string> entryLabels;
                double yMax = 0.0;
                
                if (overlayPtBins)
                {
                    for (int ipt2 = 0; ipt2 < kNPtBins; ++ipt2)
                    {
                        const PtBin& pb2 = PtBins()[ipt2];
                        const string hName = "h_ss_" + var + "_" + tag + pb2.suffix + fixedCent->suffix;
                        TH1* hSrc = GetTH1FromTopDir(aaTopSS, hName);
                        if (!hSrc) continue;
                        
                        TH1* h = CloneTH1(
                                          hSrc,
                                          TString::Format("ssQA_3x5_%s_%s_%s_%s_%s_pt%d",
                                                          H.variant.c_str(), tag.c_str(), var.c_str(),
                                                          fixedCent->folder.c_str(), trigAA.c_str(), ipt2).Data()
                                          );
                        if (!h) continue;
                        
                        EnsureSumw2(h);
                        h->SetTitle("");
                        h->SetLineColor(overlayColors[ipt2 % nOverlayColors]);
                        h->SetMarkerColor(overlayColors[ipt2 % nOverlayColors]);
                        h->SetMarkerStyle(20);
                        h->SetFillStyle(0);
                        h->SetLineWidth(2);
                        h->SetMarkerSize(0.90);
                        
                        for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                            yMax = std::max(yMax, (double)(h->GetBinContent(ib) + h->GetBinError(ib)));
                        
                        hOverlays.push_back(h);
                        entryLabels.push_back(TString::Format("%d-%d GeV", pb2.lo, pb2.hi).Data());
                        keepAlive.push_back(h);
                    }
                }
                else
                {
                    for (std::size_t ic2 = 0; ic2 < centBins.size(); ++ic2)
                    {
                        const auto& cb2 = centBins[ic2];
                        const string hName = "h_ss_" + var + "_" + tag + fixedPt->suffix + cb2.suffix;
                        TH1* hSrc = GetTH1FromTopDir(aaTopSS, hName);
                        if (!hSrc) continue;
                        
                        TH1* h = CloneNormalizeStyle(
                                                     hSrc,
                                                     TString::Format("ssQA_3x5_%s_%s_%s_%s_%s_cent%zu",
                                                                     H.variant.c_str(), tag.c_str(), var.c_str(),
                                                                     fixedPt->folder.c_str(), trigAA.c_str(), ic2).Data(),
                                                     overlayColors[(int)ic2 % nOverlayColors],
                                                     20
                                                     );
                        if (!h) continue;
                        
                        h->SetFillStyle(0);
                        h->SetLineWidth(2);
                        h->SetMarkerSize(0.90);
                        
                        for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                            yMax = std::max(yMax, (double)(h->GetBinContent(ib) + h->GetBinError(ib)));
                        
                        hOverlays.push_back(h);
                        entryLabels.push_back(TString::Format("%d-%d%%", cb2.lo, cb2.hi).Data());
                        keepAlive.push_back(h);
                    }
                }
                
                if (hOverlays.empty())
                {
                    TLatex tMiss;
                    tMiss.SetNDC(true);
                    tMiss.SetTextFont(42);
                    tMiss.SetTextAlign(22);
                    tMiss.SetTextSize(0.075);
                    tMiss.DrawLatex(0.50, 0.55, "MISSING");
                    continue;
                }
                
                anyPad = true;
                
                TH1* hFrame = hOverlays[0];
                hFrame->GetXaxis()->SetTitle(vlabel.c_str());
                hFrame->GetYaxis()->SetTitle(overlayPtBins ? "Counts" : "Unit Normalized");
                hFrame->GetYaxis()->SetTitleOffset(1.45);
                hFrame->GetYaxis()->SetTitleSize(0.050);
                hFrame->GetYaxis()->SetLabelSize(0.040);
                hFrame->SetMinimum(0.0);
                hFrame->SetMaximum((yMax > 0.0) ? (yMax * 1.10) : 1.0);
                hFrame->Draw("E1");
                for (std::size_t ih = 1; ih < hOverlays.size(); ++ih) hOverlays[ih]->Draw("E1 same");
                
                if (ivar == 0)
                {
                    TLegend* leg = new TLegend(
                                               overlayPtBins ? 0.16 : 0.18,
                                               0.52,
                                               overlayPtBins ? 0.62 : 0.56,
                                               0.88
                                               );
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(42);
                    leg->SetTextSize(overlayPtBins ? 0.024 : 0.030);
                    if (overlayPtBins) leg->SetNColumns(2);
                    
                    for (std::size_t ih = 0; ih < hOverlays.size(); ++ih)
                        leg->AddEntry(hOverlays[ih], entryLabels[ih].c_str(), "ep");
                    leg->Draw();
                    keepLegs.push_back(leg);
                }
                
                {
                    TLatex th;
                    th.SetNDC(true);
                    th.SetTextFont(42);
                    th.SetTextAlign(22);
                    th.SetTextSize(0.046);
                    th.DrawLatex(0.50, 0.91,
                                 TString::Format("%s, %s, %s",
                                                 vlabel.c_str(), TagLabel(tag).c_str(), H.label.c_str()).Data());
                }
                
                {
                    TLatex tf;
                    tf.SetNDC(true);
                    tf.SetTextFont(42);
                    tf.SetTextAlign(22);
                    tf.SetTextSize(0.036);
                    if (overlayPtBins)
                    {
                        tf.DrawLatex(0.50, 0.84,
                                     TString::Format("AuAu %d-%d%%", fixedCent->lo, fixedCent->hi).Data());
                    }
                    else
                    {
                        tf.DrawLatex(0.50, 0.84,
                                     TString::Format("p_{T}^{#gamma}: %d-%d GeV", fixedPt->lo, fixedPt->hi).Data());
                    }
                }
            }
        }
        
        if (anyPad)
        {
            SaveCanvas(cTbl, JoinPath(outDir, outName));
        }
        
        for (TLegend* leg : keepLegs) delete leg;
        for (TH1* h : keepAlive) delete h;
    };
    
    {
        const int nOverlayColorsPerPt = (int)(sizeof(overlayColors) / sizeof(overlayColors[0]));
        
        for (std::size_t ic = 0; ic < centBins.size(); ++ic)
        {
            const auto& cb = centBins[ic];
            const string centOutDir = JoinPath(perPtBinOverlayBase, cb.folder);
            EnsureDir(centOutDir);
            
            for (const auto& tag : centOverlayTags)
            {
                const string tagDir = JoinPath(centOutDir, TagFolder(tag));
                EnsureDir(tagDir);
                
                DrawSSOverlayTable3x5(
                                      tagDir,
                                      TString::Format("table3x5_SS_%s_overlayByPt.png", tag.c_str()).Data(),
                                      tag,
                                      true,
                                      &cb,
                                      nullptr
                                      );
                
                for (std::size_t vidx : ssTableVariantIdx)
                {
                    if (vidx >= handles.size()) continue;
                    auto& H = handles[vidx];
                    if (!H.file) continue;
                    
                    TDirectory* aaTopSS = H.file->GetDirectory(trigAA.c_str());
                    if (!aaTopSS) continue;
                    
                    for (int ivar = 0; ivar < (int)ssVars.size(); ++ivar)
                    {
                        const std::string& var = ssVars[ivar].var;
                        const std::string& vlabel = ssVars[ivar].label;
                        
                        const string varDir = JoinPath(centOutDir, var);
                        const string varTagDir = JoinPath(varDir, TagFolder(tag));
                        EnsureDir(varDir);
                        EnsureDir(varTagDir);
                        
                        vector<TH1*> hOverlays;
                        vector<string> ptLabels;
                        double yMax = 0.0;
                        
                        for (int ipt2 = 0; ipt2 < kNPtBins; ++ipt2)
                        {
                            const PtBin& pb2 = PtBins()[ipt2];
                            const string hName = "h_ss_" + var + "_" + tag + pb2.suffix + cb.suffix;
                            TH1* hSrc = GetTH1FromTopDir(aaTopSS, hName);
                            if (!hSrc) continue;
                            
                            TH1* h = CloneTH1(
                                              hSrc,
                                              TString::Format("ssQA_ptOverlay_%s_%s_%s_%s_%s_pt%d",
                                                              H.variant.c_str(), tag.c_str(), var.c_str(),
                                                              cb.folder.c_str(), trigAA.c_str(), ipt2).Data()
                                              );
                            if (!h) continue;
                            
                            EnsureSumw2(h);
                            h->SetTitle("");
                            h->SetLineColor(overlayColors[ipt2 % nOverlayColorsPerPt]);
                            h->SetMarkerColor(overlayColors[ipt2 % nOverlayColorsPerPt]);
                            h->SetMarkerStyle(20);
                            h->SetFillStyle(0);
                            h->SetLineWidth(2);
                            h->SetMarkerSize(0.90);
                            
                            for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                                yMax = std::max(yMax, (double)(h->GetBinContent(ib) + h->GetBinError(ib)));
                            
                            hOverlays.push_back(h);
                            ptLabels.push_back(TString::Format("%d-%d GeV", pb2.lo, pb2.hi).Data());
                        }
                        
                        if (hOverlays.empty()) continue;
                        
                        TCanvas cVarPt(
                                       TString::Format("c_ssQA_ptOv_%s_%s_%s_%s_%s",
                                                       var.c_str(), H.variant.c_str(), tag.c_str(),
                                                       cb.folder.c_str(), trigAA.c_str()).Data(),
                                       "c_ssQA_ptOv_single", 900, 700
                                       );
                        ApplyCanvasMargins1D(cVarPt);
                        cVarPt.cd();
                        
                        TH1* hFrame = hOverlays[0];
                        hFrame->GetXaxis()->SetTitle(vlabel.c_str());
                        hFrame->GetYaxis()->SetTitle("Counts");
                        hFrame->GetYaxis()->SetTitleOffset(1.15);
                        hFrame->SetMinimum(0.0);
                        hFrame->SetMaximum((yMax > 0.0) ? (yMax * 1.25) : 1.0);
                        hFrame->Draw("E1");
                        for (std::size_t ih = 1; ih < hOverlays.size(); ++ih)
                            hOverlays[ih]->Draw("E1 same");
                        hOverlays[0]->Draw("E1 same");
                        
                        TLegend legVarPt(0.2, 0.78, 0.65, 0.90);
                        legVarPt.SetBorderSize(0);
                        legVarPt.SetFillStyle(0);
                        legVarPt.SetTextFont(42);
                        legVarPt.SetTextSize(0.028);
                        legVarPt.SetNColumns(3);
                        for (std::size_t ih = 0; ih < hOverlays.size(); ++ih)
                            legVarPt.AddEntry(hOverlays[ih], ptLabels[ih].c_str(), "ep");
                        legVarPt.Draw();
                        
                        TLatex tVarPt;
                        tVarPt.SetNDC(true);
                        tVarPt.SetTextFont(42);
                        tVarPt.SetTextAlign(23);
                        tVarPt.SetTextSize(0.042);
                        tVarPt.DrawLatex(0.50, 0.955,
                                         TString::Format("%s p_{T} overlay, %s, %s",
                                                         vlabel.c_str(), TagLabel(tag).c_str(), H.label.c_str()).Data());
                        
                        TLatex tCentPt;
                        tCentPt.SetNDC(true);
                        tCentPt.SetTextFont(42);
                        tCentPt.SetTextAlign(33);
                        tCentPt.SetTextSize(0.042);
                        tCentPt.DrawLatex(0.93, 0.90,
                                          TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                        
                        DrawSSOverlayCutsAndText(var, tag, true, 0.0, 0.22, 0.735, 0.028);
                        
                        SaveCanvas(cVarPt, JoinPath(varTagDir,
                                                    TString::Format("ptOverlay_%s_%s_%s.png",
                                                                    var.c_str(), tag.c_str(), H.variant.c_str()).Data()));
                        
                        for (TH1* h : hOverlays) delete h;
                    }
                }
            }
        }
    }
    
    for (int ipt = 0; ipt < kNPtBins; ++ipt)
    {
        const PtBin& pb = PtBins()[ipt];
        const string ptOutDir = JoinPath(perCentralityOverlayBase, pb.folder);
        EnsureDir(ptOutDir);
        
        for (const auto& tag : centOverlayTags)
        {
            const string tagDir = JoinPath(ptOutDir, TagFolder(tag));
            EnsureDir(tagDir);
            
            DrawSSOverlayTable3x5(
                                  tagDir,
                                  TString::Format("table3x5_SS_%s_overlayByCent.png", tag.c_str()).Data(),
                                  tag,
                                  false,
                                  nullptr,
                                  &pb
                                  );
        }
    }
    
    // Per-variant 1x5 centrality-overlay tables (one PNG per variant per tag per pT bin)
    for (int ipt = 0; ipt < kNPtBins; ++ipt)
    {
        const PtBin& pb = PtBins()[ipt];
        const string ptOutDir = JoinPath(perCentralityOverlayBase, pb.folder);
        EnsureDir(ptOutDir);
        
        const int nOverlayColors = (int)(sizeof(overlayColors) / sizeof(overlayColors[0]));
        
        for (const auto& tag : centOverlayTags)
        {
            for (std::size_t vidx : ssTableVariantIdx)
            {
                if (vidx >= handles.size()) continue;
                auto& H = handles[vidx];
                if (!H.file) continue;
                
                TDirectory* aaTopSS = H.file->GetDirectory(trigAA.c_str());
                if (!aaTopSS) continue;
                
                TCanvas c1x5(
                             TString::Format("c_ssQA_1x5_centOv_%s_%s_%s_%s",
                                             H.variant.c_str(), tag.c_str(), pb.folder.c_str(), trigAA.c_str()).Data(),
                             "c_ssQA_1x5_centOv", 2600, 750
                             );
                c1x5.Divide(5, 1, 0.001, 0.001);
                
                vector<TH1*> keepAlive1x5;
                keepAlive1x5.reserve(ssVars.size() * centBins.size());
                
                vector<TLegend*> keepLegs1x5;
                keepLegs1x5.reserve(1);
                
                bool anyPad1x5 = false;
                
                struct VarOverlayInfo {
                    int ivar;
                    vector<TH1*> hists;
                    vector<string> labels;
                    double yMax;
                    TH1* hPP = nullptr;
                };
                vector<VarOverlayInfo> varInfos;
                varInfos.reserve(ssVars.size());
                
                for (int ivar = 0; ivar < (int)ssVars.size(); ++ivar)
                {
                    c1x5.cd(ivar + 1);
                    gPad->SetLeftMargin(0.14);
                    gPad->SetRightMargin(0.05);
                    gPad->SetBottomMargin(0.14);
                    gPad->SetTopMargin(0.18);
                    gPad->SetLogy(false);
                    
                    const std::string& var = ssVars[ivar].var;
                    const std::string& vlabel = ssVars[ivar].label;
                    
                    vector<TH1*> hOverlays;
                    vector<string> entryLabels;
                    double yMax = 0.0;
                    
                    for (std::size_t ic2 = 0; ic2 < centBins.size(); ++ic2)
                    {
                        const auto& cb2 = centBins[ic2];
                        const string hName = "h_ss_" + var + "_" + tag + pb.suffix + cb2.suffix;
                        TH1* hSrc = GetTH1FromTopDir(aaTopSS, hName);
                        if (!hSrc) continue;
                        
                        TH1* h = CloneNormalizeStyle(
                                                     hSrc,
                                                     TString::Format("ssQA_1x5_%s_%s_%s_%s_%s_cent%zu",
                                                                     H.variant.c_str(), tag.c_str(), var.c_str(),
                                                                     pb.folder.c_str(), trigAA.c_str(), ic2).Data(),
                                                     overlayColors[(int)ic2 % nOverlayColors],
                                                     20
                                                     );
                        if (!h) continue;
                        
                        h->SetFillStyle(0);
                        h->SetLineWidth(2);
                        h->SetMarkerSize(0.90);
                        
                        for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                            yMax = std::max(yMax, (double)(h->GetBinContent(ib) + h->GetBinError(ib)));
                        
                        hOverlays.push_back(h);
                        entryLabels.push_back(TString::Format("%d-%d%%", cb2.lo, cb2.hi).Data());
                        keepAlive1x5.push_back(h);
                    }
                    
                    if (hOverlays.empty())
                    {
                        TLatex tMiss;
                        tMiss.SetNDC(true);
                        tMiss.SetTextFont(42);
                        tMiss.SetTextAlign(22);
                        tMiss.SetTextSize(0.075);
                        tMiss.DrawLatex(0.50, 0.55, "MISSING");
                        continue;
                    }
                    
                    anyPad1x5 = true;
                    
                    // -- pp reference (open red circles) --
                    TH1* hPPov = nullptr;
                    if (ppTop)
                    {
                        const string hPPName = "h_ss_" + var + "_" + tag + pb.suffix;
                        TH1* hPPsrc = GetTH1FromTopDir(ppTop, hPPName);
                        if (hPPsrc)
                        {
                            hPPov = CloneNormalizeStyle(
                                                        hPPsrc,
                                                        TString::Format("ssQA_1x5_pp_%s_%s_%s_%s_%s",
                                                                        H.variant.c_str(), tag.c_str(), var.c_str(),
                                                                        pb.folder.c_str(), trigAA.c_str()).Data(),
                                                        kRed + 1, 24
                                                        );
                            if (hPPov)
                            {
                                hPPov->SetLineWidth(2);
                                hPPov->SetLineStyle(1);
                                hPPov->SetFillStyle(0);
                                hPPov->SetMarkerStyle(24);
                                hPPov->SetMarkerSize(1.00);
                                hPPov->SetMarkerColor(kRed + 1);
                                hPPov->SetLineColor(kRed + 1);
                                
                                for (int ib = 1; ib <= hPPov->GetNbinsX(); ++ib)
                                    yMax = std::max(yMax, (double)(hPPov->GetBinContent(ib) + hPPov->GetBinError(ib)));
                                
                                keepAlive1x5.push_back(hPPov);
                            }
                        }
                    }
                    
                    TH1* hFrame = hOverlays[0];
                    hFrame->GetXaxis()->SetTitle(vlabel.c_str());
                    hFrame->GetYaxis()->SetTitle("Unit Normalized");
                    hFrame->GetYaxis()->SetTitleOffset(1.45);
                    hFrame->GetYaxis()->SetTitleSize(0.050);
                    hFrame->GetYaxis()->SetLabelSize(0.040);
                    hFrame->SetMinimum(0.0);
                    hFrame->SetMaximum((yMax > 0.0) ? (yMax * 1.10) : 1.0);
                    hFrame->Draw("E1");
                    if (hPPov) hPPov->Draw("E1 same");
                    for (std::size_t ih = 1; ih < hOverlays.size(); ++ih) hOverlays[ih]->Draw("E1 same");
                    hOverlays[0]->Draw("E1 same");
                    
                    if (ivar == 0)
                    {
                        TLegend* leg = new TLegend(0.12, 0.75, 0.65, 0.92);
                        leg->SetBorderSize(0);
                        leg->SetFillStyle(0);
                        leg->SetTextFont(42);
                        leg->SetTextSize(0.037);
                        leg->SetNColumns(3);
                        
                        if (hPPov) leg->AddEntry(hPPov, "pp", "ep");
                        for (std::size_t ih = 0; ih < hOverlays.size(); ++ih)
                            leg->AddEntry(hOverlays[ih], entryLabels[ih].c_str(), "ep");
                        leg->Draw();
                        keepLegs1x5.push_back(leg);
                    }
                    
                    {
                        TLatex th;
                        th.SetNDC(true);
                        th.SetTextFont(42);
                        th.SetTextAlign(22);
                        th.SetTextSize(0.046);
                        th.DrawLatex(0.50, 0.91,
                                     TString::Format("%s, %s, %s",
                                                     vlabel.c_str(), TagLabel(tag).c_str(), H.label.c_str()).Data());
                    }
                    
                    {
                        TLatex tf;
                        tf.SetNDC(true);
                        tf.SetTextFont(42);
                        tf.SetTextAlign(33);
                        tf.SetTextSize(0.065);
                        tf.DrawLatex(0.93, 0.90,
                                     TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data());
                    }
                    
                    if (!hOverlays.empty())
                    {
                        VarOverlayInfo vi;
                        vi.ivar = ivar;
                        vi.hists = hOverlays;
                        vi.labels = entryLabels;
                        vi.yMax = yMax;
                        vi.hPP = hPPov;
                        varInfos.push_back(vi);
                    }
                }
                
                const string tagDir1x5 = JoinPath(ptOutDir, TagFolder(tag));
                EnsureDir(tagDir1x5);
                
                if (anyPad1x5)
                {
                    SaveCanvas(c1x5, JoinPath(tagDir1x5,
                                              TString::Format("table1x5_SS_%s_centOverlay_%s.png", tag.c_str(), H.variant.c_str()).Data()));
                }
                
                // --- Individual per-variable centrality overlay PNGs ---
                for (const auto& vi : varInfos)
                {
                    const std::string& var   = ssVars[vi.ivar].var;
                    const std::string& vlabel = ssVars[vi.ivar].label;
                    
                    const string varDir = JoinPath(ptOutDir, var);
                    const string varTagDir = JoinPath(varDir, TagFolder(tag));
                    EnsureDir(varTagDir);
                    
                    TCanvas cVar(
                                 TString::Format("c_ssQA_centOv_%s_%s_%s_%s_%s",
                                                 var.c_str(), H.variant.c_str(), tag.c_str(),
                                                 pb.folder.c_str(), trigAA.c_str()).Data(),
                                 "c_ssQA_centOv_single", 900, 700
                                 );
                    ApplyCanvasMargins1D(cVar);
                    cVar.cd();
                    
                    TH1* hFrame = vi.hists[0];
                    hFrame->GetXaxis()->SetTitle(vlabel.c_str());
                    hFrame->GetYaxis()->SetTitle("Unit Normalized");
                    hFrame->GetYaxis()->SetTitleOffset(1.15);
                    hFrame->SetMinimum(0.0);
                    hFrame->SetMaximum((vi.yMax > 0.0) ? (vi.yMax * 1.25) : 1.0);
                    hFrame->Draw("E1");
                    if (vi.hPP) vi.hPP->Draw("E1 same");
                    for (std::size_t ih = 1; ih < vi.hists.size(); ++ih)
                        vi.hists[ih]->Draw("E1 same");
                    vi.hists[0]->Draw("E1 same");
                    
                    const bool useSpecialE11Legend = (pb.folder == "pT_10_12" && var == "e11e33");
                    const bool useSpecialE32Legend =
                    (pb.folder == "pT_10_12" && var == "e32e35" &&
                     (tag == "pre" || tag == "tight" || tag == "nonTight"));
                    
                    TLegend legVar(
                                   useSpecialE32Legend ? 0.08 : (useSpecialE11Legend ? 0.17 : 0.56),
                                   useSpecialE32Legend ? 0.79 : (useSpecialE11Legend ? 0.78 : 0.58),
                                   useSpecialE32Legend ? 0.66 : (useSpecialE11Legend ? 0.69 : 0.92),
                                   useSpecialE32Legend ? 0.91 : (useSpecialE11Legend ? 0.90 : 0.88)
                                   );
                    legVar.SetBorderSize(0);
                    legVar.SetFillStyle(0);
                    legVar.SetTextFont(42);
                    legVar.SetTextSize((useSpecialE11Legend || useSpecialE32Legend) ? 0.028 : 0.032);
                    if (useSpecialE11Legend || useSpecialE32Legend) legVar.SetNColumns(3);
                    if (vi.hPP) legVar.AddEntry(vi.hPP, "pp", "ep");
                    for (std::size_t ih = 0; ih < vi.hists.size(); ++ih)
                        legVar.AddEntry(vi.hists[ih], vi.labels[ih].c_str(), "ep");
                    legVar.Draw();
                    
                    TLatex tVar;
                    tVar.SetNDC(true);
                    tVar.SetTextFont(42);
                    tVar.SetTextAlign(23);
                    tVar.SetTextSize(0.042);
                    tVar.DrawLatex(0.50, 0.955,
                                   TString::Format("%s centrality overlay, %s, %s",
                                                   vlabel.c_str(), TagLabel(tag).c_str(), H.label.c_str()).Data());
                    
                    TLatex tPt;
                    tPt.SetNDC(true);
                    tPt.SetTextFont(42);
                    tPt.SetTextAlign(33);
                    tPt.SetTextSize(0.042);
                    tPt.DrawLatex(0.93, 0.90,
                                  TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data());
                    
                    const double ptCenterForCuts = 0.5 * (pb.lo + pb.hi);
                    const double cutTextX = useSpecialE32Legend ? 0.08 : 0.16;
                    const double cutTextY = useSpecialE32Legend ? 0.73 : 0.84;
                    const double cutTextSize = useSpecialE32Legend ? 0.028 : 0.026;
                    DrawSSOverlayCutsAndText(var, tag, false, ptCenterForCuts, cutTextX, cutTextY, cutTextSize);
                    
                    SaveCanvas(cVar, JoinPath(varTagDir,
                                              TString::Format("centOverlay_%s_%s_%s.png",
                                                              var.c_str(), tag.c_str(), H.variant.c_str()).Data()));
                    
                    const string overlay2CentDir = JoinPath(varTagDir, "overlay2centralities");
                    EnsureDir(overlay2CentDir);
                    
                    vector<TH1*> reducedHists;
                    vector<string> reducedLabels;
                    double reducedYMax = 0.0;
                    
                    for (std::size_t ih = 0; ih < vi.hists.size(); ++ih)
                    {
                        if (vi.labels[ih] != "0-10%" && vi.labels[ih] != "60-80%") continue;
                        
                        reducedHists.push_back(vi.hists[ih]);
                        reducedLabels.push_back(vi.labels[ih]);
                        
                        TH1* hTmp = vi.hists[ih];
                        if (!hTmp) continue;
                        for (int ib = 1; ib <= hTmp->GetNbinsX(); ++ib)
                            reducedYMax = std::max(reducedYMax, (double)(hTmp->GetBinContent(ib) + hTmp->GetBinError(ib)));
                    }
                    
                    if (vi.hPP)
                    {
                        for (int ib = 1; ib <= vi.hPP->GetNbinsX(); ++ib)
                            reducedYMax = std::max(reducedYMax, (double)(vi.hPP->GetBinContent(ib) + vi.hPP->GetBinError(ib)));
                    }
                    
                    if (!reducedHists.empty())
                    {
                        TCanvas cVar2(
                                      TString::Format("c_ssQA_centOv2_%s_%s_%s_%s_%s",
                                                      var.c_str(), H.variant.c_str(), tag.c_str(),
                                                      pb.folder.c_str(), trigAA.c_str()).Data(),
                                      "c_ssQA_centOv_single_2cent", 900, 700
                                      );
                        ApplyCanvasMargins1D(cVar2);
                        cVar2.cd();
                        
                        const bool zoomOverlay2CentE32E35 =
                        (pb.folder == "pT_10_12" && var == "e32e35" && tag == "inclusive");
                        const bool useSpecialE32Legend2 =
                        (pb.folder == "pT_10_12" && var == "e32e35" && tag == "inclusive");
                        
                        TH1* hFrame2 = reducedHists[0];
                        hFrame2->GetXaxis()->SetTitle(vlabel.c_str());
                        hFrame2->GetYaxis()->SetTitle("Unit Normalized");
                        hFrame2->GetYaxis()->SetTitleOffset(1.15);
                        hFrame2->SetMinimum(0.0);
                        hFrame2->SetMaximum((reducedYMax > 0.0) ? (reducedYMax * 1.25) : 1.0);
                        if (zoomOverlay2CentE32E35) hFrame2->GetXaxis()->SetRangeUser(0.8, 1.0);
                        hFrame2->Draw("E1");
                        if (vi.hPP) vi.hPP->Draw("E1 same");
                        for (std::size_t ih = 1; ih < reducedHists.size(); ++ih)
                            reducedHists[ih]->Draw("E1 same");
                        reducedHists[0]->Draw("E1 same");
                        
                        const bool useSpecialE11Legend2 = (pb.folder == "pT_10_12" && var == "e11e33");
                        const bool useSpecialEt1Legend2 =
                        (pb.folder == "pT_10_12" && var == "et1" &&
                         (tag == "pre" || tag == "tight" || tag == "nonTight"));
                        
                        TLegend legVar2(
                                        useSpecialEt1Legend2 ? 0.08 : (useSpecialE32Legend2 ? 0.14 : (useSpecialE11Legend2 ? 0.17 : 0.56)),
                                        useSpecialEt1Legend2 ? 0.79 : (useSpecialE32Legend2 ? 0.82 : (useSpecialE11Legend2 ? 0.82 : 0.58)),
                                        useSpecialEt1Legend2 ? 0.67 : (useSpecialE32Legend2 ? 0.62 : (useSpecialE11Legend2 ? 0.69 : 0.92)),
                                        useSpecialEt1Legend2 ? 0.91 : (useSpecialE32Legend2 ? 0.90 : (useSpecialE11Legend2 ? 0.90 : 0.88))
                                        );
                        legVar2.SetBorderSize(0);
                        legVar2.SetFillStyle(0);
                        legVar2.SetTextFont(42);
                        legVar2.SetTextSize((useSpecialE11Legend2 || useSpecialE32Legend2 || useSpecialEt1Legend2) ? 0.028 : 0.032);
                        if (useSpecialE11Legend2 || useSpecialE32Legend2 || useSpecialEt1Legend2) legVar2.SetNColumns(3);
                        if (vi.hPP) legVar2.AddEntry(vi.hPP, "pp", "ep");
                        for (std::size_t ih = 0; ih < reducedHists.size(); ++ih)
                            legVar2.AddEntry(reducedHists[ih], reducedLabels[ih].c_str(), "ep");
                        legVar2.Draw();
                        
                        TLatex tVar2;
                        tVar2.SetNDC(true);
                        tVar2.SetTextFont(42);
                        tVar2.SetTextAlign(23);
                        tVar2.SetTextSize(0.042);
                        tVar2.DrawLatex(0.50, 0.955,
                                        TString::Format("%s centrality overlay, %s, %s",
                                                        vlabel.c_str(), TagLabel(tag).c_str(), H.label.c_str()).Data());
                        
                        TLatex tPt2;
                        tPt2.SetNDC(true);
                        tPt2.SetTextFont(42);
                        tPt2.SetTextAlign(33);
                        tPt2.SetTextSize(0.042);
                        tPt2.DrawLatex(0.93, 0.90,
                                       TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data());
                        
                        const double ptCenterForCuts2 = 0.5 * (pb.lo + pb.hi);
                        const double cutTextX2 = useSpecialEt1Legend2 ? 0.08 : 0.16;
                        const double cutTextY2 = useSpecialEt1Legend2 ? 0.73 : 0.84;
                        const double cutTextSize2 = useSpecialEt1Legend2 ? 0.028 : 0.026;
                        DrawSSOverlayCutsAndText(var, tag, false, ptCenterForCuts2, cutTextX2, cutTextY2, cutTextSize2);
                        
                        SaveCanvas(cVar2, JoinPath(overlay2CentDir,
                                                   TString::Format("centOverlay_%s_%s_%s.png",
                                                                   var.c_str(), tag.c_str(), H.variant.c_str()).Data()));
                    }
                }
                
                // --- Embedded Signal/Bkg MC overlays per SS var ---
                for (const auto& embSrc : embeddedSSSources)
                {
                    if (!embSrc.topDir) continue;
                    
                    for (const auto& vi : varInfos)
                    {
                        const std::string& var   = ssVars[vi.ivar].var;
                        const std::string& vlabel = ssVars[vi.ivar].label;
                        
                        const string varDir = JoinPath(ptOutDir, var);
                        const string varTagDir = JoinPath(varDir, TagFolder(tag));
                        const string embDir = JoinPath(varTagDir, embSrc.folderName);
                        EnsureDir(embDir);
                        
                        // Fetch signal / bkg MC from embedded sim
                        const string hSigName = "h_ss_" + var + "_" + tag + "_sig" + pb.suffix;
                        const string hBkgName = "h_ss_" + var + "_" + tag + "_bkg" + pb.suffix;
                        
                        TH1* rawSig = GetTH1FromTopDir(embSrc.topDir, hSigName);
                        TH1* rawBkg = GetTH1FromTopDir(embSrc.topDir, hBkgName);
                        
                        TH1* hSig = nullptr;
                        TH1* hBkg = nullptr;
                        vector<TH1*> embKeep;
                        
                        double embYMax = vi.yMax;
                        
                        if (rawSig)
                        {
                            hSig = CloneNormalizeStyle(
                                rawSig,
                                TString::Format("ssQA_centOv_embSig_%s_%s_%s_%s_%s_%s",
                                    embSrc.folderName.c_str(), var.c_str(), tag.c_str(),
                                    H.variant.c_str(), pb.folder.c_str(), trigAA.c_str()).Data(),
                                kPink + 7, 1);
                            if (hSig)
                            {
                                hSig->SetLineWidth(2);
                                hSig->SetLineColor(kPink + 7);
                                hSig->SetMarkerStyle(1);
                                hSig->SetMarkerSize(0.0);
                                hSig->SetFillStyle(0);
                                for (int ib = 1; ib <= hSig->GetNbinsX(); ++ib)
                                    embYMax = std::max(embYMax, (double)(hSig->GetBinContent(ib) + hSig->GetBinError(ib)));
                                embKeep.push_back(hSig);
                            }
                        }
                        
                        if (rawBkg)
                        {
                            hBkg = CloneNormalizeStyle(
                                rawBkg,
                                TString::Format("ssQA_centOv_embBkg_%s_%s_%s_%s_%s_%s",
                                    embSrc.folderName.c_str(), var.c_str(), tag.c_str(),
                                    H.variant.c_str(), pb.folder.c_str(), trigAA.c_str()).Data(),
                                kBlue + 1, 1);
                            if (hBkg)
                            {
                                hBkg->SetLineWidth(2);
                                hBkg->SetLineColor(kBlue + 1);
                                hBkg->SetMarkerStyle(1);
                                hBkg->SetMarkerSize(0.0);
                                hBkg->SetFillStyle(0);
                                for (int ib = 1; ib <= hBkg->GetNbinsX(); ++ib)
                                    embYMax = std::max(embYMax, (double)(hBkg->GetBinContent(ib) + hBkg->GetBinError(ib)));
                                embKeep.push_back(hBkg);
                            }
                        }
                        
                        // --- Full centrality overlay + embedded MC ---
                        {
                            TCanvas cEmb(
                                TString::Format("c_ssQA_centOv_emb_%s_%s_%s_%s_%s_%s",
                                    embSrc.folderName.c_str(), var.c_str(), H.variant.c_str(),
                                    tag.c_str(), pb.folder.c_str(), trigAA.c_str()).Data(),
                                "c_ssQA_centOv_emb", 900, 700);
                            ApplyCanvasMargins1D(cEmb);
                            cEmb.cd();
                            
                            TH1* hFrame = vi.hists[0];
                            hFrame->GetXaxis()->SetTitle(vlabel.c_str());
                            hFrame->GetYaxis()->SetTitle("Unit Normalized");
                            hFrame->GetYaxis()->SetTitleOffset(1.15);
                            hFrame->SetMinimum(0.0);
                            hFrame->SetMaximum((embYMax > 0.0) ? (embYMax * 1.25) : 1.0);
                            hFrame->Draw("E1");
                            
                            if (hSig) hSig->Draw("HIST same");
                            if (hBkg) hBkg->Draw("HIST same");
                            if (vi.hPP) vi.hPP->Draw("E1 same");
                            for (std::size_t ih = 1; ih < vi.hists.size(); ++ih)
                                vi.hists[ih]->Draw("E1 same");
                            vi.hists[0]->Draw("E1 same");
                            
                            const bool useSpecialE11Leg = (pb.folder == "pT_10_12" && var == "e11e33");
                            const bool useSpecialE32Leg =
                                (pb.folder == "pT_10_12" && var == "e32e35" &&
                                 (tag == "pre" || tag == "tight" || tag == "nonTight"));
                            
                            TLegend legEmb(
                                useSpecialE32Leg ? 0.08 : (useSpecialE11Leg ? 0.17 : 0.56),
                                useSpecialE32Leg ? 0.74 : (useSpecialE11Leg ? 0.72 : 0.50),
                                useSpecialE32Leg ? 0.66 : (useSpecialE11Leg ? 0.69 : 0.92),
                                useSpecialE32Leg ? 0.91 : (useSpecialE11Leg ? 0.90 : 0.88));
                            legEmb.SetBorderSize(0);
                            legEmb.SetFillStyle(0);
                            legEmb.SetTextFont(42);
                            legEmb.SetTextSize((useSpecialE11Leg || useSpecialE32Leg) ? 0.026 : 0.028);
                            if (useSpecialE11Leg || useSpecialE32Leg) legEmb.SetNColumns(3);
                            if (vi.hPP) legEmb.AddEntry(vi.hPP, "pp", "ep");
                            for (std::size_t ih = 0; ih < vi.hists.size(); ++ih)
                                legEmb.AddEntry(vi.hists[ih], vi.labels[ih].c_str(), "ep");
                            if (hSig) legEmb.AddEntry(hSig, "Signal MC", "l");
                            if (hBkg) legEmb.AddEntry(hBkg, "Background MC", "l");
                            legEmb.Draw();
                            
                            TLatex tEmb;
                            tEmb.SetNDC(true);
                            tEmb.SetTextFont(42);
                            tEmb.SetTextAlign(23);
                            tEmb.SetTextSize(0.038);
                            tEmb.DrawLatex(0.50, 0.955,
                                TString::Format("%s centrality overlay, %s, %s (%s)",
                                    vlabel.c_str(), TagLabel(tag).c_str(), H.label.c_str(),
                                    embSrc.label.c_str()).Data());
                            
                            TLatex tPtEmb;
                            tPtEmb.SetNDC(true);
                            tPtEmb.SetTextFont(42);
                            tPtEmb.SetTextAlign(33);
                            tPtEmb.SetTextSize(0.042);
                            tPtEmb.DrawLatex(0.93, 0.90,
                                TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data());
                            
                            const double ptCE = 0.5 * (pb.lo + pb.hi);
                            const double cutXE = useSpecialE32Leg ? 0.08 : 0.16;
                            const double cutYE = useSpecialE32Leg ? 0.68 : 0.84;
                            const double cutSE = useSpecialE32Leg ? 0.028 : 0.026;
                            DrawSSOverlayCutsAndText(var, tag, false, ptCE, cutXE, cutYE, cutSE);
                            
                            SaveCanvas(cEmb, JoinPath(embDir,
                                TString::Format("centOverlay_%s_%s_%s.png",
                                    var.c_str(), tag.c_str(), H.variant.c_str()).Data()));
                        }
                        
                        // --- overlay2centralities (0-10% + 60-80%) + embedded MC ---
                        {
                            vector<TH1*> redH;
                            vector<string> redL;
                            double redYMax = 0.0;
                            
                            for (std::size_t ih = 0; ih < vi.hists.size(); ++ih)
                            {
                                if (vi.labels[ih] != "0-10%" && vi.labels[ih] != "60-80%") continue;
                                redH.push_back(vi.hists[ih]);
                                redL.push_back(vi.labels[ih]);
                                TH1* hT = vi.hists[ih];
                                for (int ib = 1; ib <= hT->GetNbinsX(); ++ib)
                                    redYMax = std::max(redYMax, (double)(hT->GetBinContent(ib) + hT->GetBinError(ib)));
                            }
                            if (vi.hPP)
                                for (int ib = 1; ib <= vi.hPP->GetNbinsX(); ++ib)
                                    redYMax = std::max(redYMax, (double)(vi.hPP->GetBinContent(ib) + vi.hPP->GetBinError(ib)));
                            if (hSig)
                                for (int ib = 1; ib <= hSig->GetNbinsX(); ++ib)
                                    redYMax = std::max(redYMax, (double)(hSig->GetBinContent(ib) + hSig->GetBinError(ib)));
                            if (hBkg)
                                for (int ib = 1; ib <= hBkg->GetNbinsX(); ++ib)
                                    redYMax = std::max(redYMax, (double)(hBkg->GetBinContent(ib) + hBkg->GetBinError(ib)));
                            
                            if (!redH.empty())
                            {
                                const string ov2Dir = JoinPath(embDir, "overlay2centralities");
                                EnsureDir(ov2Dir);
                                
                                TCanvas cEmb2(
                                    TString::Format("c_ssQA_centOv2_emb_%s_%s_%s_%s_%s_%s",
                                        embSrc.folderName.c_str(), var.c_str(), H.variant.c_str(),
                                        tag.c_str(), pb.folder.c_str(), trigAA.c_str()).Data(),
                                    "c_ssQA_centOv2_emb", 900, 700);
                                ApplyCanvasMargins1D(cEmb2);
                                cEmb2.cd();
                                
                                TH1* hFrame2 = redH[0];
                                hFrame2->GetXaxis()->SetTitle(vlabel.c_str());
                                hFrame2->GetYaxis()->SetTitle("Unit Normalized");
                                hFrame2->GetYaxis()->SetTitleOffset(1.15);
                                hFrame2->SetMinimum(0.0);
                                hFrame2->SetMaximum((redYMax > 0.0) ? (redYMax * 1.25) : 1.0);
                                hFrame2->Draw("E1");
                                
                                if (hSig) hSig->Draw("HIST same");
                                if (hBkg) hBkg->Draw("HIST same");
                                if (vi.hPP) vi.hPP->Draw("E1 same");
                                for (std::size_t ih = 1; ih < redH.size(); ++ih)
                                    redH[ih]->Draw("E1 same");
                                redH[0]->Draw("E1 same");
                                
                                TLegend legEmb2(0.56, 0.50, 0.92, 0.88);
                                legEmb2.SetBorderSize(0);
                                legEmb2.SetFillStyle(0);
                                legEmb2.SetTextFont(42);
                                legEmb2.SetTextSize(0.028);
                                if (vi.hPP) legEmb2.AddEntry(vi.hPP, "pp", "ep");
                                for (std::size_t ih = 0; ih < redH.size(); ++ih)
                                    legEmb2.AddEntry(redH[ih], redL[ih].c_str(), "ep");
                                if (hSig) legEmb2.AddEntry(hSig, "Signal MC", "l");
                                if (hBkg) legEmb2.AddEntry(hBkg, "Background MC", "l");
                                legEmb2.Draw();
                                
                                TLatex tE2;
                                tE2.SetNDC(true);
                                tE2.SetTextFont(42);
                                tE2.SetTextAlign(23);
                                tE2.SetTextSize(0.038);
                                tE2.DrawLatex(0.50, 0.955,
                                    TString::Format("%s centrality overlay, %s, %s (%s)",
                                        vlabel.c_str(), TagLabel(tag).c_str(), H.label.c_str(),
                                        embSrc.label.c_str()).Data());
                                
                                TLatex tPE2;
                                tPE2.SetNDC(true);
                                tPE2.SetTextFont(42);
                                tPE2.SetTextAlign(33);
                                tPE2.SetTextSize(0.042);
                                tPE2.DrawLatex(0.93, 0.90,
                                    TString::Format("p_{T}^{#gamma}: %d-%d GeV", pb.lo, pb.hi).Data());
                                
                                const double ptCE2 = 0.5 * (pb.lo + pb.hi);
                                DrawSSOverlayCutsAndText(var, tag, false, ptCE2, 0.16, 0.84, 0.026);
                                
                                SaveCanvas(cEmb2, JoinPath(ov2Dir,
                                    TString::Format("centOverlay_%s_%s_%s.png",
                                        var.c_str(), tag.c_str(), H.variant.c_str()).Data()));
                            }
                        }
                        
                        for (TH1* h : embKeep) delete h;
                    }
                }
                
                for (TLegend* leg : keepLegs1x5) delete leg;
                for (TH1* h : keepAlive1x5) delete h;
            }
        }
    }
    
    // Close embedded SS overlay files
    for (auto& src : embeddedSSSources)
    {
        src.topDir = nullptr;
        if (src.file) { src.file->Close(); delete src.file; src.file = nullptr; }
    }
}

{
  const string passFracBase = JoinPath(ssQADir, "passFractions");
  EnsureDir(passFracBase);

  struct PFCategoryDef
  {
    string key;
    string histTag;
  };

  const vector<PFCategoryDef> pfCats = {
    {"inclusive", "inclusive"},
    {"iso", "iso"},
    {"nonIso", "nonIso"},
    {"pre", "pre"},
    {"tight", "tight"},
    {"nonTight", "nonTight"}
  };

  const vector<string> pfPassDefs = {"pre", "tight"};

  struct PFSeries
  {
    string label;
    int color = kBlack;
    int marker = 20;
    vector<double> x;
    vector<double> ex;
    vector<double> y;
    vector<double> ey;
  };

  for (const auto& sv : ssVars)
  {
    const string varDir = JoinPath(passFracBase, sv.var);
    EnsureDir(varDir);

    for (const auto& passDef : pfPassDefs)
    {
      if (!PassDefValidForVar(sv.var, passDef)) continue;

      const string passDir = JoinPath(varDir, passDef);
      EnsureDir(passDir);

      for (const auto& cat : pfCats)
      {
        const string catDir = JoinPath(passDir, cat.key);
        const string perCentBase = JoinPath(catDir, "perCentrality");
        const string perPtBase   = JoinPath(catDir, "perPt");
        EnsureDir(catDir);
        EnsureDir(perCentBase);
        EnsureDir(perPtBase);

        // --------------------------------------------------
        // perCentrality/<cent>/passFraction_vsPt.png
        // --------------------------------------------------
        for (std::size_t ic = 0; ic < centBins.size(); ++ic)
        {
          const auto& cb = centBins[ic];
          const string outCentDir = JoinPath(perCentBase, cb.folder);
          EnsureDir(outCentDir);

          PFSeries sPP;
          sPP.label = "pp";
          sPP.color = kRed + 1;
          sPP.marker = 24;

            PFSeries sNoSub;
            sNoSub.label = TString::Format("AuAu (%d-%d%%) No UE sub", cb.lo, cb.hi).Data();
            sNoSub.color = ssColors[0];
            sNoSub.marker = ssMarkers[0];

            PFSeries sVarA;
            sVarA.label = TString::Format("AuAu (%d-%d%%) Variant A", cb.lo, cb.hi).Data();
            sVarA.color = ssColors[2];
            sVarA.marker = ssMarkers[2];

            PFSeries sVarB;
            sVarB.label = TString::Format("AuAu (%d-%d%%) Variant B", cb.lo, cb.hi).Data();
            sVarB.color = ssColors[3];
            sVarB.marker = ssMarkers[3];

          for (int ipt = 0; ipt < kNPtBins; ++ipt)
          {
            const PtBin& pb = PtBins()[ipt];
            const double ptCenter = 0.5 * (pb.lo + pb.hi);
            const double ptHalfW  = 0.5 * (pb.hi - pb.lo);

            const string hPPName = "h_ss_" + sv.var + "_" + cat.histTag + pb.suffix;
            TH1* hPPsrc = GetTH1FromTopDir(ppTop, hPPName);
            if (hPPsrc)
            {
              double frac = 0.0, efrac = 0.0;
              if (PassFractionForHist(hPPsrc, sv.var, passDef, ptCenter, frac, efrac))
              {
                sPP.x.push_back(ptCenter);
                sPP.ex.push_back(ptHalfW);
                sPP.y.push_back(frac);
                sPP.ey.push_back(efrac);
              }
            }

            auto FillVariantPoint =
              [&](std::size_t idx, PFSeries& S)
            {
              if (idx >= handles.size()) return;
              auto& H = handles[idx];
              if (!H.file) return;

              TDirectory* aaTopSS = H.file->GetDirectory(trigAA.c_str());
              if (!aaTopSS) return;

              const string hAAName = "h_ss_" + sv.var + "_" + cat.histTag + pb.suffix + cb.suffix;
              TH1* hAAsrc = GetTH1FromTopDir(aaTopSS, hAAName);
              if (!hAAsrc) return;

              double frac = 0.0, efrac = 0.0;
              if (PassFractionForHist(hAAsrc, sv.var, passDef, ptCenter, frac, efrac))
              {
                S.x.push_back(ptCenter);
                S.ex.push_back(ptHalfW);
                S.y.push_back(frac);
                S.ey.push_back(efrac);
              }
            };

            FillVariantPoint(0, sNoSub);
            FillVariantPoint(2, sVarA);
            FillVariantPoint(3, sVarB);
          }

          const bool haveAny =
            !sPP.x.empty() || !sNoSub.x.empty() || !sVarA.x.empty() || !sVarB.x.empty();

          if (!haveAny) continue;

          TCanvas cPF(
            TString::Format("c_pf_vsPt_%s_%s_%s_%s_%s",
              sv.var.c_str(), passDef.c_str(), cat.key.c_str(), cb.folder.c_str(), trigAA.c_str()).Data(),
            "c_pf_vsPt", 900, 700
          );
          ApplyCanvasMargins1D(cPF);
          cPF.cd();

          TH1F hFrame(
            TString::Format("hFrame_pf_vsPt_%s_%s_%s_%s",
              sv.var.c_str(), passDef.c_str(), cat.key.c_str(), cb.folder.c_str()).Data(),
            "", 100, kPtEdges.front(), kPtEdges.back()
          );
          hFrame.SetDirectory(nullptr);
          hFrame.SetStats(0);
          hFrame.SetMinimum(0.0);
          hFrame.SetMaximum(1.05);
          hFrame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
          hFrame.GetYaxis()->SetTitle("Pass fraction");
          hFrame.GetXaxis()->SetTitleSize(0.055);
          hFrame.GetYaxis()->SetTitleSize(0.055);
          hFrame.GetXaxis()->SetLabelSize(0.045);
          hFrame.GetYaxis()->SetLabelSize(0.045);
          hFrame.GetYaxis()->SetTitleOffset(1.15);
          hFrame.Draw();

          vector<TGraphErrors*> keepGraphs;

          auto DrawSeries =
            [&](const PFSeries& S) -> TGraphErrors*
            {
              if (S.x.empty()) return nullptr;
              TGraphErrors* g = new TGraphErrors(
                (int)S.x.size(),
                const_cast<double*>(&S.x[0]),
                const_cast<double*>(&S.y[0]),
                const_cast<double*>(&S.ex[0]),
                const_cast<double*>(&S.ey[0])
              );
              g->SetLineWidth(2);
              g->SetLineColor(S.color);
              g->SetMarkerColor(S.color);
              g->SetMarkerStyle(S.marker);
              g->SetMarkerSize(1.2);
              g->Draw("PE1 SAME");
              return g;
            };

          TGraphErrors* gPP    = DrawSeries(sPP);
          TGraphErrors* gNoSub = DrawSeries(sNoSub);
          TGraphErrors* gVarA  = DrawSeries(sVarA);
          TGraphErrors* gVarB  = DrawSeries(sVarB);

          if (gPP)    keepGraphs.push_back(gPP);
          if (gNoSub) keepGraphs.push_back(gNoSub);
          if (gVarA)  keepGraphs.push_back(gVarA);
          if (gVarB)  keepGraphs.push_back(gVarB);

            TLegend leg(0.44, 0.66, 0.80, 0.88);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.SetTextFont(42);
            leg.SetTextSize(0.032);
            if (gPP)    leg.AddEntry(gPP,    sPP.label.c_str(), "ep");
            if (gNoSub) leg.AddEntry(gNoSub, sNoSub.label.c_str(), "ep");
            if (gVarA)  leg.AddEntry(gVarA,  sVarA.label.c_str(), "ep");
            if (gVarB)  leg.AddEntry(gVarB,  sVarB.label.c_str(), "ep");
            leg.Draw();

          TLatex tTitle;
          tTitle.SetNDC(true);
          tTitle.SetTextFont(42);
          tTitle.SetTextAlign(23);
          tTitle.SetTextSize(0.042);
          tTitle.DrawLatex(0.50, 0.965,
            TString::Format("%s %s pass fraction vs p_{T}^{#gamma}, %s, %d-%d%%",
              sv.label.c_str(), passDef.c_str(), cat.key.c_str(), cb.lo, cb.hi).Data());

          SaveCanvas(cPF, JoinPath(outCentDir, "passFraction_vsPt.png"));

          for (auto* g : keepGraphs) delete g;
        }

        // --------------------------------------------------
        // perPt/<pTbin>/passFraction_vsCent.png
        // --------------------------------------------------
        for (int ipt = 0; ipt < kNPtBins; ++ipt)
        {
          const PtBin& pb = PtBins()[ipt];
          const string outPtDir = JoinPath(perPtBase, pb.folder);
          EnsureDir(outPtDir);

          PFSeries sPP;
          sPP.label = "pp";
          sPP.color = kRed + 1;
          sPP.marker = 24;

          PFSeries sNoSub;
          sNoSub.label = "No UE sub";
          sNoSub.color = ssColors[0];
          sNoSub.marker = ssMarkers[0];

          PFSeries sVarA;
          sVarA.label = "Variant A";
          sVarA.color = ssColors[2];
          sVarA.marker = ssMarkers[2];

          PFSeries sVarB;
          sVarB.label = "Variant B";
          sVarB.color = ssColors[3];
          sVarB.marker = ssMarkers[3];

          const double ptCenter = 0.5 * (pb.lo + pb.hi);

          for (std::size_t ic = 0; ic < centBins.size(); ++ic)
          {
            const auto& cb = centBins[ic];
            const double centCenter = 0.5 * (cb.lo + cb.hi);
            const double centHalfW  = 0.5 * (cb.hi - cb.lo);

            const string hPPName = "h_ss_" + sv.var + "_" + cat.histTag + pb.suffix;
            TH1* hPPsrc = GetTH1FromTopDir(ppTop, hPPName);
            if (hPPsrc)
            {
              double frac = 0.0, efrac = 0.0;
              if (PassFractionForHist(hPPsrc, sv.var, passDef, ptCenter, frac, efrac))
              {
                sPP.x.push_back(centCenter);
                sPP.ex.push_back(centHalfW);
                sPP.y.push_back(frac);
                sPP.ey.push_back(efrac);
              }
            }

            auto FillVariantPoint =
              [&](std::size_t idx, PFSeries& S)
            {
              if (idx >= handles.size()) return;
              auto& H = handles[idx];
              if (!H.file) return;

              TDirectory* aaTopSS = H.file->GetDirectory(trigAA.c_str());
              if (!aaTopSS) return;

              const string hAAName = "h_ss_" + sv.var + "_" + cat.histTag + pb.suffix + cb.suffix;
              TH1* hAAsrc = GetTH1FromTopDir(aaTopSS, hAAName);
              if (!hAAsrc) return;

              double frac = 0.0, efrac = 0.0;
              if (PassFractionForHist(hAAsrc, sv.var, passDef, ptCenter, frac, efrac))
              {
                S.x.push_back(centCenter);
                S.ex.push_back(centHalfW);
                S.y.push_back(frac);
                S.ey.push_back(efrac);
              }
            };

            FillVariantPoint(0, sNoSub);
            FillVariantPoint(2, sVarA);
            FillVariantPoint(3, sVarB);
          }

          const bool haveAny =
            !sPP.x.empty() || !sNoSub.x.empty() || !sVarA.x.empty() || !sVarB.x.empty();

          if (!haveAny) continue;

          TCanvas cPF(
            TString::Format("c_pf_vsCent_%s_%s_%s_%s_%s",
              sv.var.c_str(), passDef.c_str(), cat.key.c_str(), pb.folder.c_str(), trigAA.c_str()).Data(),
            "c_pf_vsCent", 900, 700
          );
          ApplyCanvasMargins1D(cPF);
          cPF.cd();

          TH1F hFrame(
            TString::Format("hFrame_pf_vsCent_%s_%s_%s_%s",
              sv.var.c_str(), passDef.c_str(), cat.key.c_str(), pb.folder.c_str()).Data(),
            "", 100, centBins.front().lo, centBins.back().hi
          );
          hFrame.SetDirectory(nullptr);
          hFrame.SetStats(0);
          hFrame.SetMinimum(0.0);
          hFrame.SetMaximum(1.05);
          hFrame.GetXaxis()->SetTitle("Centrality [%]");
          hFrame.GetYaxis()->SetTitle("Pass fraction");
          hFrame.GetXaxis()->SetTitleSize(0.055);
          hFrame.GetYaxis()->SetTitleSize(0.055);
          hFrame.GetXaxis()->SetLabelSize(0.045);
          hFrame.GetYaxis()->SetLabelSize(0.045);
          hFrame.GetYaxis()->SetTitleOffset(1.15);
          hFrame.Draw();

          vector<TGraphErrors*> keepGraphs;

          auto DrawSeries =
            [&](const PFSeries& S) -> TGraphErrors*
            {
              if (S.x.empty()) return nullptr;
              TGraphErrors* g = new TGraphErrors(
                (int)S.x.size(),
                const_cast<double*>(&S.x[0]),
                const_cast<double*>(&S.y[0]),
                const_cast<double*>(&S.ex[0]),
                const_cast<double*>(&S.ey[0])
              );
              g->SetLineWidth(2);
              g->SetLineColor(S.color);
              g->SetMarkerColor(S.color);
              g->SetMarkerStyle(S.marker);
              g->SetMarkerSize(1.2);
              g->Draw("PE1 SAME");
              return g;
            };

          TGraphErrors* gPP    = DrawSeries(sPP);
          TGraphErrors* gNoSub = DrawSeries(sNoSub);
          TGraphErrors* gVarA  = DrawSeries(sVarA);
          TGraphErrors* gVarB  = DrawSeries(sVarB);

          if (gPP)    keepGraphs.push_back(gPP);
          if (gNoSub) keepGraphs.push_back(gNoSub);
          if (gVarA)  keepGraphs.push_back(gVarA);
          if (gVarB)  keepGraphs.push_back(gVarB);

          TLegend leg(0.56, 0.66, 0.92, 0.88);
          leg.SetBorderSize(0);
          leg.SetFillStyle(0);
          leg.SetTextFont(42);
          leg.SetTextSize(0.032);
          if (gPP)    leg.AddEntry(gPP,    sPP.label.c_str(), "ep");
          if (gNoSub) leg.AddEntry(gNoSub, sNoSub.label.c_str(), "ep");
          if (gVarA)  leg.AddEntry(gVarA,  sVarA.label.c_str(), "ep");
          if (gVarB)  leg.AddEntry(gVarB,  sVarB.label.c_str(), "ep");
          leg.Draw();

          TLatex tTitle;
          tTitle.SetNDC(true);
          tTitle.SetTextFont(42);
          tTitle.SetTextAlign(23);
          tTitle.SetTextSize(0.042);
          tTitle.DrawLatex(0.50, 0.965,
            TString::Format("%s %s pass fraction vs centrality, %s, p_{T}^{#gamma} %d-%d GeV",
              sv.label.c_str(), passDef.c_str(), cat.key.c_str(), pb.lo, pb.hi).Data());

          SaveCanvas(cPF, JoinPath(outPtDir, "passFraction_vsCent.png"));

          for (auto* g : keepGraphs) delete g;
        }
      }
    }
  }
}

if (fSimSS)
{
  fSimSS->Close();
  delete fSimSS;
  fSimSS = nullptr;
}

