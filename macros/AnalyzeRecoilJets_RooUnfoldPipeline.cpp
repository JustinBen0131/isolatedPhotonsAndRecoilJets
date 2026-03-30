#if ARJ_HAVE_ROOUNFOLD
    // =============================================================================
    // Section 5I: RooUnfold pipeline (SIM+DATA PP only)
    //   - Unfold N_{#gamma}(pT^{#gamma}) using SIM photon response
    //   - Unfold (pT^{#gamma}, x_{J}) using SIM global-bin response per radius
    //   - Produce per-photon normalized particle-level x_{J} (3x3 table: 9 analysis pT bins)
    //
    // Output base (DATA trigger):
    //   <triggerName>/unfolding/radii/<rKey>/
    //     table3x3_unfolded_perPhoton_dNdXJ.png
    //     rooUnfold_outputs.root
    //     summary_rooUnfold_pipeline.txt
    //
    // Shared photon outputs:
    //   <triggerName>/unfolding/radii/photons/
    // =============================================================================
    void RunRooUnfoldPipeline_SimAndDataPP(Dataset& dsData, Dataset& dsSim)
      {
        cout << ANSI_BOLD_CYN << "\n==============================\n"
             << "[SECTION 5I] RooUnfold pipeline (SIM+DATA PP)\n"
             << "==============================" << ANSI_RESET << "\n";

        cout << "  Compile-time: ARJ_HAVE_ROOUNFOLD=" << ARJ_HAVE_ROOUNFOLD << "\n";

        cout << "  DATA:\n"
             << "    label      = " << dsData.label << "\n"
             << "    isSim      = " << (dsData.isSim ? "true" : "false") << "\n"
             << "    trigger    = " << dsData.trigger << "\n"
             << "    topDirName = " << dsData.topDirName << "\n"
             << "    inFilePath = " << dsData.inFilePath << "\n"
             << "    outBase    = " << dsData.outBase << "\n";

        cout << "  SIM:\n"
             << "    label      = " << dsSim.label << "\n"
             << "    isSim      = " << (dsSim.isSim ? "true" : "false") << "\n"
             << "    topDirName = " << dsSim.topDirName << "\n"
             << "    inFilePath = " << dsSim.inFilePath << "\n"
             << "    outBase    = " << dsSim.outBase << "\n";

        if (!dsData.file || !dsData.topDir || !dsSim.file || !dsSim.topDir)
        {
          cout << ANSI_BOLD_RED
               << "[ERROR] RooUnfold pipeline cannot run because dataset file/topDir is null.\n"
               << "        DATA: file=" << dsData.file << " topDir=" << dsData.topDir << "\n"
               << "        SIM : file=" << dsSim.file  << " topDir=" << dsSim.topDir  << "\n"
               << ANSI_RESET << "\n";
          return;
        }

        cout << "  DATA topDir path: " << dsData.topDir->GetPath() << "\n";
        cout << "  SIM  topDir path: " << dsSim.topDir->GetPath()  << "\n";

        if (gSystem)
        {
          const int rc = gSystem->Load("libRooUnfold");
          cout << "  gSystem->Load(\"libRooUnfold\") rc=" << rc << "\n";
          if (rc < 0)
          {
            cout << "  gSystem->GetDynamicPath() = " << gSystem->GetDynamicPath() << "\n";
            cout << ANSI_BOLD_YEL
                 << "[WARN] gSystem->Load(\"libRooUnfold\") failed. Headers were found (compiled), but runtime RooUnfold symbols may be missing."
                 << ANSI_RESET << "\n";
          }
      }

      if (dsData.isSim || !dsSim.isSim)
      {
            cout << ANSI_BOLD_YEL << "[WARN] RunRooUnfoldPipeline_SimAndDataPP called with unexpected dataset types (dsData.isSim="
                 << dsData.isSim << ", dsSim.isSim=" << dsSim.isSim << "). Continuing anyway." << ANSI_RESET << "\n";
      }

        const string unfoldVariant = (gApplyPurityCorrectionForUnfolding ? "purityCorrected" : "nonPurityCorrected");
        const string combVariant   = (gApplyCombinatoricSubtractionForUnfolding ? "combinatoricSub" : "nonCombSub");
        const string simDataBase = JoinPath(dsSim.outBase, dsData.trigger);
        const bool isEmbeddedAuAu = IsEmbeddedSimSample(CurrentSimSample());

        string unfoldFolderLabel;
        if (isEmbeddedAuAu)
        {
          if (gApplyCombinatoricSubtractionForUnfolding)
          {
            unfoldFolderLabel = gApplyPurityCorrectionForUnfolding
                                  ? "combinatoricJetSubtractectedWithPurityCorrection"
                                  : "combinatoricJetSubtractectedNoPurityCorrection";
          }
          else
          {
            unfoldFolderLabel = gApplyPurityCorrectionForUnfolding
                                  ? "purityCorrected"
                                  : "nonPurityCorrected";
          }
        }
        else
        {
          unfoldFolderLabel = unfoldVariant;
        }

        const string outBase = JoinPath(simDataBase, "unfolding/" + unfoldFolderLabel + "/radii");
        int mkrc = -999;
        if (gSystem) mkrc = gSystem->mkdir(outBase.c_str(), true);
        else EnsureDir(outBase);

        cout << "  Unfolding variant: " << unfoldVariant << "\n";
        if (isEmbeddedAuAu)
        {
          cout << "  Combinatoric variant: " << combVariant
               << "  |  folder=" << unfoldFolderLabel << "\n";
        }
        cout << "  Output base: " << outBase;
        if (mkrc != -999) cout << "  (mkdir rc=" << mkrc << ")";
        cout << "\n";

        LeakageFactors lfForPurity;
        if (gApplyPurityCorrectionForUnfolding)
        {
            LoadLeakageFactorsFromSIM(dsSim, lfForPurity);
        }

        // ----------------------------------------------------------------------
        // Load ATLAS pp (HEPData ins1694678 Table 1) x_{J#gamma} points for a simple overlay
        //   CSV layout:
        //     - file contains multiple blocks; we only read the pp block after "#: Centrality,,,pp"
        //     - columns (pp block):
        //         0: xJg (center)
        //         1: xJg LOW
        //         2: xJg HIGH
        //         3: (1/N_{#gamma}) dN/dxJg in pp
        //         4: stat +
        //         6: sys +
        //   y-error used here (simple first comparison):
        //     sqrt(stat^2 + sys^2)  (symmetric)
        // ----------------------------------------------------------------------
        const string kAtlasTable1PhoPtLabel = "63.1-79.6 GeV";
        TGraphAsymmErrors* gAtlasPP = nullptr;
        string atlasCsvUsed = "";

        {
          const vector<string> csvCandidates =
          {
            "/Users/patsfan753/Desktop/ThesisAnalysis/LHC_overlayData/HEPData-ins1694678-v1-Table_1.csv",
            "LHC_overlayData/HEPData-ins1694678-v1-Table_1.csv",
            "../LHC_overlayData/HEPData-ins1694678-v1-Table_1.csv"
          };

          std::ifstream fin;
          for (const auto& p : csvCandidates)
          {
            fin.open(p);
            if (fin.good())
            {
              atlasCsvUsed = p;
              break;
            }
            fin.clear();
          }

          if (!fin.good())
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] ATLAS overlay disabled: could not open HEPData CSV (Table 1): " << csvCandidates[0] << "\n"
                 << "       (also tried: " << csvCandidates[1] << " and " << csvCandidates[2] << ")\n"
                 << ANSI_RESET << "\n";
          }
          else
          {
            vector<double> vx, vxl, vxh, vy, vstatp, vsysp;

            bool inPP = false;
            bool sawHeader = false;
            string line;

            while (std::getline(fin, line))
            {
              if (!line.empty() && line.back() == '\r') line.pop_back();

              if (line.rfind("#:", 0) == 0)
              {
                if (line.find("Centrality") != string::npos && line.find("pp") != string::npos)
                {
                  inPP = true;
                  sawHeader = false;
                }
                continue;
              }

              if (!inPP) continue;

              if (line.empty()) break;

              if (!sawHeader)
              {
                sawHeader = true;
                continue;
              }

              std::stringstream ss(line);
              vector<string> cols;
              string cell;

              while (std::getline(ss, cell, ','))
              {
                cell.erase(std::remove(cell.begin(), cell.end(), '\"'), cell.end());
                while (!cell.empty() && (cell.front() == ' ' || cell.front() == '\t')) cell.erase(cell.begin());
                while (!cell.empty() && (cell.back()  == ' ' || cell.back()  == '\t')) cell.pop_back();
                cols.push_back(cell);
              }

              if (cols.size() < 8) continue;

              auto toD = [](const string& s, double& out) -> bool
              {
                try
                {
                  size_t idx = 0;
                  out = std::stod(s, &idx);
                  return (idx > 0);
                }
                catch (...)
                {
                  return false;
                }
              };

              double x = 0.0, xlo = 0.0, xhi = 0.0, y = 0.0, statp = 0.0, sysp = 0.0;

              if (!toD(cols[0], x))   continue;
              if (!toD(cols[1], xlo)) continue;
              if (!toD(cols[2], xhi)) continue;
              if (!toD(cols[3], y))   continue;

              if (!toD(cols[4], statp)) statp = 0.0;
              if (!toD(cols[6], sysp))  sysp  = 0.0;

              vx.push_back(x);
              vxl.push_back(xlo);
              vxh.push_back(xhi);
              vy.push_back(y);
              vstatp.push_back(statp);
              vsysp.push_back(sysp);
            }

            if (!vx.empty())
            {
              gAtlasPP = new TGraphAsymmErrors((int)vx.size());
              gAtlasPP->SetName("gAtlas_pp_Table1");
              gAtlasPP->SetTitle("");

              for (int i = 0; i < (int)vx.size(); ++i)
              {
                const double exl = vx[i] - vxl[i];
                const double exh = vxh[i] - vx[i];
                const double ey  = std::sqrt(vstatp[i]*vstatp[i] + vsysp[i]*vsysp[i]);

                gAtlasPP->SetPoint(i, vx[i], vy[i]);
                gAtlasPP->SetPointError(i, exl, exh, ey, ey);
              }

              gAtlasPP->SetMarkerStyle(25);
              gAtlasPP->SetMarkerSize(1.00);
              gAtlasPP->SetMarkerColor(2);
              gAtlasPP->SetLineColor(2);
              gAtlasPP->SetLineWidth(2);

              cout << "  ATLAS overlay loaded (pp, Table 1): " << atlasCsvUsed
                   << "  [pT^gamma=" << kAtlasTable1PhoPtLabel << "]  N=" << gAtlasPP->GetN() << "\n";
            }
            else
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] ATLAS overlay disabled: parsed zero pp points from CSV: " << atlasCsvUsed
                   << ANSI_RESET << "\n";
            }
          }
        }

      // ----------------------------------------------------------------------
      // Helper: transpose any TH2 (including variable binning) so that:
      //   hOut(x=oldY, y=oldX) with all bin contents/errors preserved
      // ----------------------------------------------------------------------
      auto TransposeTH2 = [](const TH2* hIn, const string& newName, const string& newTitle) -> TH2D*
      {
        if (!hIn) return nullptr;

        const TAxis* ax = hIn->GetXaxis();
        const TAxis* ay = hIn->GetYaxis();

        const int nx = ax->GetNbins();
        const int ny = ay->GetNbins();

        const bool xVar = (ay->GetXbins() && ay->GetXbins()->GetSize() > 0);
        const bool yVar = (ax->GetXbins() && ax->GetXbins()->GetSize() > 0);

        TH2D* hOut = nullptr;

        if (xVar && yVar)
        {
          hOut = new TH2D(newName.c_str(), newTitle.c_str(),
                          ny, ay->GetXbins()->GetArray(),
                          nx, ax->GetXbins()->GetArray());
        }
        else if (xVar && !yVar)
        {
          hOut = new TH2D(newName.c_str(), newTitle.c_str(),
                          ny, ay->GetXbins()->GetArray(),
                          nx, ax->GetXmin(), ax->GetXmax());
        }
        else if (!xVar && yVar)
        {
          hOut = new TH2D(newName.c_str(), newTitle.c_str(),
                          ny, ay->GetXmin(), ay->GetXmax(),
                          nx, ax->GetXbins()->GetArray());
        }
        else
        {
          hOut = new TH2D(newName.c_str(), newTitle.c_str(),
                          ny, ay->GetXmin(), ay->GetXmax(),
                          nx, ax->GetXmin(), ax->GetXmax());
        }

        hOut->SetDirectory(nullptr);
        hOut->Sumw2();

        // Include underflow/overflow so nothing is silently dropped.
        for (int ix = 0; ix <= nx + 1; ++ix)
        {
          for (int iy = 0; iy <= ny + 1; ++iy)
          {
            hOut->SetBinContent(iy, ix, hIn->GetBinContent(ix, iy));
            hOut->SetBinError  (iy, ix, hIn->GetBinError  (ix, iy));
          }
        }

        return hOut;
      };

      // ----------------------------------------------------------------------
      // Helper: flatten TH2(pT,xJ) into TH1(globalBin) using ROOT's internal
      //         global-bin indexing:
      //           g = (nx+2)*iy + ix, with ix,iy including under/overflow.
      //
      // This is the SAME g used when filling h2_unfoldResponse_* in RecoilJets.cc:
      //   gReco   = h2Reco ->FindBin(pT_reco, xJ_reco);
      //   gTruth  = h2Truth->FindBin(pT_true, xJ_true);
      // ----------------------------------------------------------------------
      auto FlattenTH2ToGlobal = [](const TH2* h2, const string& newName) -> TH1D*
      {
        if (!h2) return nullptr;

        const int nx = h2->GetNbinsX();
        const int ny = h2->GetNbinsY();
        const int nGlob = (nx + 2) * (ny + 2);

        TH1D* h1 = new TH1D(newName.c_str(), "", nGlob, -0.5, nGlob - 0.5);
        h1->SetDirectory(nullptr);
        h1->Sumw2();

        for (int ix = 0; ix <= nx + 1; ++ix)
        {
          for (int iy = 0; iy <= ny + 1; ++iy)
          {
            const int g = h2->GetBin(ix, iy);     // 0 .. nGlob-1
            const int b = g + 1;                  // TH1 bin number (1-based)
            if (b < 1 || b > nGlob) continue;

            h1->SetBinContent(b, h2->GetBinContent(ix, iy));
            h1->SetBinError  (b, h2->GetBinError  (ix, iy));
          }
        }

        return h1;
      };

      // ----------------------------------------------------------------------
      // Helper: un-flatten TH1(globalTruth) back into TH2(truth pT,xJ) using a
      //         template TH2 for binning.
      // ----------------------------------------------------------------------
        auto UnflattenGlobalToTH2 = [](const TH1* hGlob, const TH2* tmpl, const string& newName) -> TH2*
        {
          if (!hGlob || !tmpl) return nullptr;

          TH2* h2 = CloneTH2(tmpl, newName);
          if (!h2) return nullptr;

          h2->Reset("ICES");
          h2->SetDirectory(nullptr);

          const int nx = h2->GetNbinsX();
          const int ny = h2->GetNbinsY();
          const int nGlob = (nx + 2) * (ny + 2);

          if (hGlob->GetNbinsX() != nGlob)
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] UnflattenGlobalToTH2: template nGlob=" << nGlob
                 << " but hGlob nbins=" << hGlob->GetNbinsX()
                 << ". Proceeding (bins beyond overlap will be ignored)."
                 << ANSI_RESET << "\n";
          }

          const int nCopy = std::min(nGlob, hGlob->GetNbinsX());

          for (int ix = 0; ix <= nx + 1; ++ix)
          {
            for (int iy = 0; iy <= ny + 1; ++iy)
            {
              const int g = h2->GetBin(ix, iy); // 0..nGlob-1
              const int b = g + 1;
              if (b < 1 || b > nCopy) continue;

              h2->SetBinContent(ix, iy, hGlob->GetBinContent(b));
              h2->SetBinError  (ix, iy, hGlob->GetBinError  (b));
            }
          }

          return h2;
        };

        // ----------------------------------------------------------------------
        // Helpers for ATLAS-style purity correction before unfolding:
        //   - photon normalization input: use S_A(pT^{#gamma}) from ABCD
        //   - measured xJ input: H_A - (N_{bkg}^{A} / N_C) * H_C  in each pT^{#gamma} bin
        //
        // The response matrix remains the nominal signal response (truth signal -> reco A),
        // which is the correct thing to use once the fake-photon recoil template has been
        // subtracted from the measured reco distribution before unfolding.
        // ----------------------------------------------------------------------
        auto ComputeABCDSignalCounts =
          [&](int iPt,
              double& A,
              double& B,
              double& C,
              double& D,
              double& SA,
              double& eSA,
              double& nbkgA)->void
        {
          const PtBin& b = PtBins()[iPt];

          A = Read1BinCount(dsData, "h_isIsolated_isTight"   + b.suffix);
          B = Read1BinCount(dsData, "h_notIsolated_isTight"  + b.suffix);
          C = Read1BinCount(dsData, "h_isIsolated_notTight"  + b.suffix);
          D = Read1BinCount(dsData, "h_notIsolated_notTight" + b.suffix);

          SA = 0.0;
          if (A > 0.0 && D > 0.0)
          {
            SA = A - B * (C / D);
            if (SA < 0.0) SA = 0.0;
          }

          if (gApplyPurityCorrectionForUnfolding && lfForPurity.available)
          {
            double SAcorr = 0.0;
            if (SolveLeakageCorrectedSA(A, B, C, D,
                                        lfForPurity.fB[iPt],
                                        lfForPurity.fC[iPt],
                                        lfForPurity.fD[iPt],
                                        SAcorr))
            {
              SA = std::min(std::max(SAcorr, 0.0), A);
            }
          }

          double varSA = 0.0;
          if (A > 0.0) varSA += 1.0 * 1.0 * A;
          if (D > 0.0)
          {
            const double dSdB = -(C / D);
            const double dSdC = -(B / D);
            const double dSdD =  (B * C) / (D * D);

            if (B > 0.0) varSA += dSdB * dSdB * B;
            if (C > 0.0) varSA += dSdC * dSdC * C;
            if (D > 0.0) varSA += dSdD * dSdD * D;
          }

          eSA = (varSA > 0.0) ? std::sqrt(varSA) : 0.0;
          nbkgA = std::max(0.0, A - SA);
        };

        auto ApplyPurityCorrectionToRecoPhotonHist =
          [&](TH1* h)->bool
        {
          if (!h) return false;

          h->Reset("ICES");
          h->SetDirectory(nullptr);
          EnsureSumw2(h);

          const auto& recoBins = UnfoldRecoPtBins();
          const int nRecoBins = (int)recoBins.size();

          for (int i = 0; i < nRecoBins; ++i)
          {
            double A = 0.0, B = 0.0, C = 0.0, D = 0.0;
            double SA = 0.0, eSA = 0.0, nbkgA = 0.0;

            const PtBin& b = recoBins[i];

            if (b.lo == 8  && b.hi == 10)
            {
              ComputeABCDSignalCounts(0, A, B, C, D, SA, eSA, nbkgA);
            }
            else if (b.lo == 35 && b.hi == 40)
            {
              ComputeABCDSignalCounts(kNPtBins - 1, A, B, C, D, SA, eSA, nbkgA);
            }
            else
            {
              int iCanon = -1;
              for (int j = 0; j < kNPtBins; ++j)
              {
                const PtBin& pb = PtBins()[j];
                if (pb.lo == b.lo && pb.hi == b.hi)
                {
                  iCanon = j;
                  break;
                }
              }

              if (iCanon < 0) continue;
              ComputeABCDSignalCounts(iCanon, A, B, C, D, SA, eSA, nbkgA);
            }

            const double cen = 0.5 * (b.lo + b.hi);
            const int ib = h->GetXaxis()->FindBin(cen);
            if (ib < 1 || ib > h->GetNbinsX()) continue;

            h->SetBinContent(ib, SA);
            h->SetBinError  (ib, eSA);
          }

          return true;
        };

        auto ApplyPurityCorrectionToRecoXJHist =
          [&](TH2* h, TH2* hSideC)->bool
        {
          if (!h || !hSideC) return false;

          TH2* hAorig = CloneTH2(h, TString::Format("%s_Aorig_forPurity", h->GetName()).Data());
          if (!hAorig) return false;
          hAorig->SetDirectory(nullptr);
          EnsureSumw2(hAorig);

          h->Reset("ICES");
          h->SetDirectory(nullptr);
          EnsureSumw2(h);

          const int ny = h->GetYaxis()->GetNbins();
          const auto& recoBins = UnfoldRecoPtBins();
          const int nRecoBins = (int)recoBins.size();

          for (int i = 0; i < nRecoBins; ++i)
          {
            double A = 0.0, B = 0.0, C = 0.0, D = 0.0;
            double SA = 0.0, eSA = 0.0, nbkgA = 0.0;

            const PtBin& b = recoBins[i];

            if (b.lo == 8  && b.hi == 10)
            {
              ComputeABCDSignalCounts(0, A, B, C, D, SA, eSA, nbkgA);
            }
            else if (b.lo == 35 && b.hi == 40)
            {
              ComputeABCDSignalCounts(kNPtBins - 1, A, B, C, D, SA, eSA, nbkgA);
            }
            else
            {
              int iCanon = -1;
              for (int j = 0; j < kNPtBins; ++j)
              {
                const PtBin& pb = PtBins()[j];
                if (pb.lo == b.lo && pb.hi == b.hi)
                {
                  iCanon = j;
                  break;
                }
              }

              if (iCanon < 0) continue;
              ComputeABCDSignalCounts(iCanon, A, B, C, D, SA, eSA, nbkgA);
            }

            const double scaleC = (C > 0.0) ? (nbkgA / C) : 0.0;

            const double cen = 0.5 * (b.lo + b.hi);
            const int ix = h->GetXaxis()->FindBin(cen);
            if (ix < 1 || ix > h->GetXaxis()->GetNbins()) continue;

            for (int iy = 0; iy <= ny + 1; ++iy)
            {
              const double valA = hAorig->GetBinContent(ix, iy);
              const double errA = hAorig->GetBinError  (ix, iy);

              const double valC = hSideC->GetBinContent(ix, iy);
              const double errC = hSideC->GetBinError  (ix, iy);

              double val  = valA - scaleC * valC;
              double err2 = errA * errA + scaleC * scaleC * errC * errC;

              if (!std::isfinite(val))  val  = 0.0;
              if (!std::isfinite(err2)) err2 = 0.0;

              if (val < 0.0) val = 0.0;

              h->SetBinContent(ix, iy, val);
              h->SetBinError  (ix, iy, (err2 > 0.0) ? std::sqrt(err2) : 0.0);
            }
          }

          delete hAorig;
          return true;
        };

        // ----------------------------------------------------------------------
        // (A) Photon unfolding: N_gamma(pT^gamma) at particle (truth) level
        // ----------------------------------------------------------------------

        // Shared y-axis ranges (taken from Step-B r04 closure/half-closure) to enforce identical ranges in Step-A photon plots
        static double gYmin_r04_closure     = std::numeric_limits<double>::quiet_NaN();
        static double gYmax_r04_closure     = std::numeric_limits<double>::quiet_NaN();
        static double gYmin_r04_halfClosure = std::numeric_limits<double>::quiet_NaN();
        static double gYmax_r04_halfClosure = std::numeric_limits<double>::quiet_NaN();

        const string phoDir = JoinPath(outBase, "photons");
        EnsureDir(phoDir);
        const string phoYieldDir = JoinPath(phoDir, "yields");
        EnsureDir(phoYieldDir);

        auto MakePhotonFigure30Like =
          [&](TH1* hRecoDataForFit,
              TH1* hRecoSimForFit,
              const string& outYieldDir,
              const string& plotTag,
              const string& recoDefinitionLabel)->void
        {
          if (!hRecoDataForFit || !hRecoSimForFit) return;

          const double fitMin = 8.0;
          const double fitMax = 35.0;

          if (hRecoDataForFit->GetNbinsX() != hRecoSimForFit->GetNbinsX())
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] Figure-30-like photon fit skipped (" << plotTag
                 << "): reco DATA/SIM histograms do not have the same number of bins."
                 << ANSI_RESET << "\n";
            return;
          }

          struct BinPoint
          {
            double lo = 0.0;
            double hi = 0.0;
            double x = 0.0;
            double ex = 0.0;
            double data = 0.0;
            double edata = 0.0;
            double sim = 0.0;
            double esim = 0.0;
          };

          vector<BinPoint> bins;
          bins.reserve((size_t) hRecoDataForFit->GetNbinsX());

          double sumData = 0.0;
          double sumSim  = 0.0;

          for (int ib = 1; ib <= hRecoDataForFit->GetNbinsX(); ++ib)
          {
            const double dLo = hRecoDataForFit->GetXaxis()->GetBinLowEdge(ib);
            const double dHi = hRecoDataForFit->GetXaxis()->GetBinUpEdge(ib);
            const double sLo = hRecoSimForFit ->GetXaxis()->GetBinLowEdge(ib);
            const double sHi = hRecoSimForFit ->GetXaxis()->GetBinUpEdge(ib);

            if (std::fabs(dLo - sLo) > 1e-9 || std::fabs(dHi - sHi) > 1e-9)
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Figure-30-like photon fit skipped (" << plotTag
                   << "): reco DATA/SIM histogram bin edges do not match."
                   << ANSI_RESET << "\n";
              return;
            }

            if (dLo < fitMin || dHi > fitMax) continue;

            BinPoint p;
            p.lo = dLo;
            p.hi = dHi;
            p.x  = 0.5 * (dLo + dHi);
            p.ex = 0.5 * (dHi - dLo);

            p.data  = hRecoDataForFit->GetBinContent(ib);
            p.edata = hRecoDataForFit->GetBinError  (ib);
            p.sim   = hRecoSimForFit ->GetBinContent(ib);
            p.esim  = hRecoSimForFit ->GetBinError  (ib);

            if (p.edata <= 0.0 && p.data > 0.0) p.edata = std::sqrt(p.data);
            if (p.esim  <= 0.0 && p.sim  > 0.0) p.esim  = std::sqrt(p.sim);

            sumData += std::max(0.0, p.data);
            sumSim  += std::max(0.0, p.sim);

            bins.push_back(p);
          }

          if (bins.empty() || !(sumData > 0.0) || !(sumSim > 0.0))
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] Figure-30-like photon fit skipped (" << plotTag
                 << "): no valid 8-35 GeV reco bins with positive normalization."
                 << ANSI_RESET << "\n";
            return;
          }

          vector<double> displayEdges;
          displayEdges.reserve(bins.size() + 1);
          for (size_t i = 0; i < bins.size(); ++i) displayEdges.push_back(bins[i].lo);
          displayEdges.push_back(bins.back().hi);

          vector<double> gxData, gexData, gyData, geyData;
          vector<double> gxSim,  gexSim,  gySim,  geySim;
          vector<double> gxRatio, gexRatio, gyRatio, geyRatio;
          vector<double> gxFit,   gexFit,   gyFit,   geyFit;

          double topMax = 0.0;
          double topMinPos = std::numeric_limits<double>::max();
          double ratioMax = 1.40;

          for (size_t i = 0; i < bins.size(); ++i)
          {
            const double width = bins[i].hi - bins[i].lo;
            if (!(width > 0.0)) continue;

            const double yData = bins[i].data / (sumData * width);
            const double eData = bins[i].edata / (sumData * width);
            const double ySim  = bins[i].sim  / (sumSim  * width);
            const double eSim  = bins[i].esim / (sumSim  * width);

            gxData.push_back(bins[i].x);
            gexData.push_back(bins[i].ex);
            gyData.push_back(yData);
            geyData.push_back(eData);

            gxSim.push_back(bins[i].x);
            gexSim.push_back(bins[i].ex);
            gySim.push_back(ySim);
            geySim.push_back(eSim);

            if (std::isfinite(yData) && yData > 0.0) topMinPos = std::min(topMinPos, yData);
            if (std::isfinite(ySim)  && ySim  > 0.0) topMinPos = std::min(topMinPos, ySim);

            if (std::isfinite(yData + eData)) topMax = std::max(topMax, yData + eData);
            if (std::isfinite(ySim  + eSim )) topMax = std::max(topMax, ySim  + eSim );

            double ratio = 0.0;
            double eratio = 0.0;

            if (ySim > 0.0)
            {
              ratio = yData / ySim;

              if (yData > 0.0)
              {
                const double relData = (eData > 0.0) ? (eData / yData) : 0.0;
                const double relSim  = (eSim  > 0.0) ? (eSim  / ySim ) : 0.0;
                eratio = ratio * std::sqrt(relData * relData + relSim * relSim);
              }
              else
              {
                eratio = (eData > 0.0) ? (eData / ySim) : 0.0;
              }
            }

            gxRatio.push_back(bins[i].x);
            gexRatio.push_back(bins[i].ex);
            gyRatio.push_back(ratio);
            geyRatio.push_back(eratio);

            if (std::isfinite(ratio) && std::isfinite(eratio))
            {
              ratioMax = std::max(ratioMax, ratio + eratio);
            }

            if (ratio > 0.0 && eratio > 0.0 && std::isfinite(ratio) && std::isfinite(eratio))
            {
              gxFit.push_back(bins[i].x);
              gexFit.push_back(bins[i].ex);
              gyFit.push_back(ratio);
              geyFit.push_back(eratio);
            }
          }

          if (!(topMinPos < std::numeric_limits<double>::max())) topMinPos = 1e-6;
          if (!(topMax > 0.0)) topMax = 1.0;

          TGraphErrors gData((int)gxData.size(), &gxData[0], &gyData[0], &gexData[0], &geyData[0]);
          TGraphErrors gSim ((int)gxSim .size(), &gxSim [0], &gySim [0], &gexSim [0], &geySim [0]);
          TGraphErrors gRatio((int)gxRatio.size(), &gxRatio[0], &gyRatio[0], &gexRatio[0], &geyRatio[0]);

          gData.SetName(TString::Format("gPhoFig30Data_%s", plotTag.c_str()).Data());
          gSim .SetName(TString::Format("gPhoFig30Sim_%s",  plotTag.c_str()).Data());
          gRatio.SetName(TString::Format("gPhoFig30Ratio_%s", plotTag.c_str()).Data());

          gData.SetMarkerStyle(20);
          gData.SetMarkerSize(1.15);
          gData.SetMarkerColor(kBlack);
          gData.SetLineColor(kBlack);
          gData.SetLineWidth(2);

          gSim.SetMarkerStyle(24);
          gSim.SetMarkerSize(1.20);
          gSim.SetMarkerColor(kMagenta + 1);
          gSim.SetLineColor(kMagenta + 1);
          gSim.SetLineWidth(2);

          gRatio.SetMarkerStyle(20);
          gRatio.SetMarkerSize(1.00);
          gRatio.SetMarkerColor(kBlack);
          gRatio.SetLineColor(kBlack);
          gRatio.SetLineWidth(2);

          TF1 fPade(
            TString::Format("fPhoFig30Pade_%s", plotTag.c_str()).Data(),
            "([0] + [1]*((x-20.0)/10.0) + [2]*((x-20.0)/10.0)^2) / (1.0 + [3]*((x-20.0)/10.0) + [4]*((x-20.0)/10.0)^2)",
            fitMin,
            fitMax
          );
          fPade.SetLineColor(kBlack);
          fPade.SetLineWidth(3);
          fPade.SetNpx(600);
          fPade.SetParNames("p0", "p1", "p2", "q1", "q2");
          fPade.SetParameters(1.0, 0.0, 0.0, 0.0, 0.0);
          fPade.SetParLimits(0, 0.0, 5.0);
          fPade.SetParLimits(1, -10.0, 10.0);
          fPade.SetParLimits(2, -10.0, 10.0);
          fPade.SetParLimits(3, -10.0, 10.0);
          fPade.SetParLimits(4, -10.0, 10.0);

          bool haveFit = false;
          int fitStatus = -999;
          double fitChi2 = 0.0;
          int fitNdf = 0;
          double fitProb = 0.0;
          TFitResultPtr fitResPtr;
          TFitResult* fitRes = nullptr;

          if (gxFit.size() >= 5)
          {
            TGraphErrors gFit((int)gxFit.size(), &gxFit[0], &gyFit[0], &gexFit[0], &geyFit[0]);
            fitResPtr = gFit.Fit(&fPade, "SQRN");
            fitStatus = (int) fitResPtr;
            fitRes = fitResPtr.Get();

            if (fitRes && fitStatus == 0)
            {
              haveFit = true;
              fitChi2 = fitRes->Chi2();
              fitNdf  = fitRes->Ndf();
              fitProb = TMath::Prob(fitChi2, fitNdf);
            }
          }

          auto EvalFrozenWeight = [&](double x)->double
          {
            double xx = x;
            if (xx < fitMin) xx = fitMin;
            if (xx > fitMax) xx = fitMax;
            return fPade.Eval(xx);
          };

          vector<double> gxPull, gexPull, gyPull, geyPull;
          double pullAbsMax = 3.0;

          if (haveFit)
          {
            for (size_t i = 0; i < bins.size(); ++i)
            {
              const double ratio = (i < gyRatio.size()) ? gyRatio[i] : 0.0;
              const double eratio = (i < geyRatio.size()) ? geyRatio[i] : 0.0;
              if (!(eratio > 0.0) || !std::isfinite(ratio) || !std::isfinite(eratio)) continue;

              const double fitVal = EvalFrozenWeight(bins[i].x);
              if (!std::isfinite(fitVal)) continue;

              const double pull = (ratio - fitVal) / eratio;

              gxPull.push_back(bins[i].x);
              gexPull.push_back(bins[i].ex);
              gyPull.push_back(pull);
              geyPull.push_back(0.0);

              if (std::isfinite(pull)) pullAbsMax = std::max(pullAbsMax, std::fabs(pull));
            }
          }

          if (pullAbsMax > 5.0) pullAbsMax = 5.0;

          TGraphErrors* gPull = nullptr;
          if (!gxPull.empty())
          {
            gPull = new TGraphErrors((int)gxPull.size(), &gxPull[0], &gyPull[0], &gexPull[0], &geyPull[0]);
            gPull->SetName(TString::Format("gPhoFig30Pull_%s", plotTag.c_str()).Data());
            gPull->SetMarkerStyle(20);
            gPull->SetMarkerSize(0.95);
            gPull->SetMarkerColor(kBlack);
            gPull->SetLineColor(kBlack);
            gPull->SetLineWidth(2);
          }

          EnsureDir(outYieldDir);

          TCanvas c(TString::Format("c_phoFig30_%s", plotTag.c_str()).Data(),
                    TString::Format("c_phoFig30_%s", plotTag.c_str()).Data(),
                    950, 1050);

          TPad* pTop = new TPad(TString::Format("pTop_phoFig30_%s", plotTag.c_str()).Data(), "",
                                0.0, 0.46, 1.0, 1.0);
          TPad* pMid = new TPad(TString::Format("pMid_phoFig30_%s", plotTag.c_str()).Data(), "",
                                0.0, 0.18, 1.0, 0.46);
          TPad* pBot = new TPad(TString::Format("pBot_phoFig30_%s", plotTag.c_str()).Data(), "",
                                0.0, 0.00, 1.0, 0.18);

          pTop->SetLeftMargin(0.14);
          pTop->SetRightMargin(0.04);
          pTop->SetTopMargin(0.06);
          pTop->SetBottomMargin(0.02);
          pTop->SetTicks(1, 1);
          pTop->SetLogy(1);

          pMid->SetLeftMargin(0.14);
          pMid->SetRightMargin(0.04);
          pMid->SetTopMargin(0.02);
          pMid->SetBottomMargin(0.02);
          pMid->SetTicks(1, 1);

          pBot->SetLeftMargin(0.14);
          pBot->SetRightMargin(0.04);
          pBot->SetTopMargin(0.02);
          pBot->SetBottomMargin(0.32);
          pBot->SetTicks(1, 1);
          pBot->SetGridy(1);

          pTop->Draw();
          pMid->Draw();
          pBot->Draw();

          pTop->cd();

          TH1D hFrameTop(TString::Format("hPhoFig30Top_%s", plotTag.c_str()).Data(),
                         "",
                         (int) displayEdges.size() - 1,
                         &displayEdges[0]);
          hFrameTop.SetDirectory(nullptr);
          hFrameTop.SetStats(0);
          hFrameTop.SetTitle("");
          hFrameTop.SetMinimum(std::max(1e-6, 0.5 * topMinPos));
          hFrameTop.SetMaximum(std::max(1.0, 1.35 * topMax));
          hFrameTop.GetXaxis()->SetTitle("");
          hFrameTop.GetXaxis()->SetLabelSize(0.0);
          hFrameTop.GetXaxis()->SetTitleSize(0.0);
          hFrameTop.GetYaxis()->SetTitle("dN/N/dE_{T}");
          hFrameTop.GetYaxis()->SetTitleSize(0.055);
          hFrameTop.GetYaxis()->SetTitleOffset(1.10);
          hFrameTop.GetYaxis()->SetLabelSize(0.045);
          hFrameTop.Draw("axis");

          gData.Draw("PZ same");
          gSim.Draw("PZ same");

          TLegend legTop(0.14, 0.14, 0.42, 0.28);
          legTop.SetBorderSize(0);
          legTop.SetFillStyle(0);
          legTop.SetTextFont(42);
          legTop.SetTextSize(0.040);
          legTop.AddEntry(&gData, "Data (purity corrected)", "pe");
          legTop.AddEntry(&gSim,  "Pythia signal", "pe");
          legTop.Draw();

          TLatex tTop;
          tTop.SetNDC(true);
          tTop.SetTextFont(42);
          tTop.SetTextAlign(33);
          tTop.SetTextSize(0.048);
          tTop.DrawLatex(0.93, 0.93, "sPHENIX Internal");
          tTop.DrawLatex(0.93, 0.86, "p+p  #sqrt{s} = 200 GeV");
          tTop.DrawLatex(0.93, 0.79, "|#eta^{#gamma}| < 0.7");
          tTop.SetTextSize(0.035);
          tTop.DrawLatex(0.93, 0.70, recoDefinitionLabel.c_str());

          pMid->cd();

          TH1D hFrameMid(TString::Format("hPhoFig30Mid_%s", plotTag.c_str()).Data(),
                         "",
                         (int) displayEdges.size() - 1,
                         &displayEdges[0]);
          hFrameMid.SetDirectory(nullptr);
          hFrameMid.SetStats(0);
          hFrameMid.SetTitle("");
          hFrameMid.SetMinimum(0.0);
          hFrameMid.SetMaximum(std::max(1.45, 1.15 * ratioMax));
          hFrameMid.GetXaxis()->SetTitle("");
          hFrameMid.GetXaxis()->SetLabelSize(0.0);
          hFrameMid.GetXaxis()->SetTitleSize(0.0);
          hFrameMid.GetYaxis()->SetTitle("Reweighting factor");
          hFrameMid.GetYaxis()->SetTitleSize(0.055);
          hFrameMid.GetYaxis()->SetTitleOffset(1.10);
          hFrameMid.GetYaxis()->SetLabelSize(0.045);
          hFrameMid.Draw("axis");

          gRatio.Draw("PZ same");
          if (haveFit) fPade.Draw("same");

          {
            TLine l1(fitMin, 1.0, fitMax, 1.0);
            l1.SetLineStyle(2);
            l1.SetLineWidth(2);
            l1.Draw("same");
          }

          TLatex tMid;
          tMid.SetNDC(true);
          tMid.SetTextFont(42);
          tMid.SetTextAlign(13);
          tMid.SetTextSize(0.038);
          tMid.DrawLatex(0.16, 0.88, "Fit: Pade [2/2], 8-35 GeV");
          tMid.DrawLatex(0.16, 0.78, "Non-positive ratio points shown, but excluded from fit");

          pBot->cd();

          TH1D hFrameBot(TString::Format("hPhoFig30Bot_%s", plotTag.c_str()).Data(),
                         "",
                         (int) displayEdges.size() - 1,
                         &displayEdges[0]);
          hFrameBot.SetDirectory(nullptr);
          hFrameBot.SetStats(0);
          hFrameBot.SetTitle("");
          hFrameBot.SetMinimum(-pullAbsMax);
          hFrameBot.SetMaximum(+pullAbsMax);
          hFrameBot.GetXaxis()->SetTitle("E_{T}^{#gamma} [GeV]");
          hFrameBot.GetXaxis()->SetTitleSize(0.15);
          hFrameBot.GetXaxis()->SetTitleOffset(0.95);
          hFrameBot.GetXaxis()->SetLabelSize(0.12);
          hFrameBot.GetYaxis()->SetTitle("Pull");
          hFrameBot.GetYaxis()->SetTitleSize(0.12);
          hFrameBot.GetYaxis()->SetTitleOffset(0.45);
          hFrameBot.GetYaxis()->SetLabelSize(0.10);
          hFrameBot.GetYaxis()->SetNdivisions(505);
          hFrameBot.Draw("axis");

          {
            TLine l0(fitMin, 0.0, fitMax, 0.0);
            l0.SetLineStyle(2);
            l0.SetLineWidth(2);
            l0.Draw("same");
          }

          if (gPull) gPull->Draw("PZ same");

          c.cd();
          c.Modified();
          c.Update();
          SaveCanvas(c, JoinPath(outYieldDir, "figure30_like_purityCorrectedData_vs_pythiaSignal.png"));

          vector<string> fitTxt;
          fitTxt.push_back("Figure-30-like photon reweighting input");
          fitTxt.push_back("Reco definition: " + recoDefinitionLabel);
          fitTxt.push_back("Displayed / fit range: 8-35 GeV");
          fitTxt.push_back("Top panel quantity: dN/(N dE_{T}) with N taken over the displayed 8-35 GeV bins");
          fitTxt.push_back("Data input: purity-corrected reco photon yield");
          fitTxt.push_back("MC input: signal reco photon yield");
          fitTxt.push_back("Fit form: Pade [2/2] in t = (E_{T}^{#gamma} - 20 GeV) / 10 GeV");
          fitTxt.push_back("w(t) = (p0 + p1 t + p2 t^{2}) / (1 + q1 t + q2 t^{2})");
          fitTxt.push_back("Evaluation rule outside fit range for future reweighting use: freeze to nearest boundary value");
          fitTxt.push_back("Non-positive ratio points are excluded from the fit but still drawn on the plot");
          fitTxt.push_back(TString::Format("N(display bins) = %d", (int) bins.size()).Data());
          fitTxt.push_back(TString::Format("Data normalization sum = %.10g", sumData).Data());
          fitTxt.push_back(TString::Format("SIM normalization sum = %.10g", sumSim).Data());
          fitTxt.push_back("");

          if (haveFit && fitRes)
          {
            fitTxt.push_back(TString::Format("Fit status = %d", fitStatus).Data());
            fitTxt.push_back(TString::Format("chi2 = %.10g", fitChi2).Data());
            fitTxt.push_back(TString::Format("ndf = %d", fitNdf).Data());
            fitTxt.push_back(TString::Format("prob = %.10g", fitProb).Data());
            fitTxt.push_back("");

            for (int ip = 0; ip < 5; ++ip)
            {
              fitTxt.push_back(TString::Format("%s = %.10g +/- %.10g",
                                               fPade.GetParName(ip),
                                               fPade.GetParameter(ip),
                                               fPade.GetParError(ip)).Data());
            }

            fitTxt.push_back("");
            fitTxt.push_back("Covariance matrix:");
            const TMatrixDSym cov = fitRes->GetCovarianceMatrix();
            for (int ir = 0; ir < cov.GetNrows(); ++ir)
            {
              std::ostringstream row;
              row << "  row " << ir << ":";
              for (int ic = 0; ic < cov.GetNcols(); ++ic)
              {
                row << " " << std::setprecision(10) << cov(ir, ic);
              }
              fitTxt.push_back(row.str());
            }
          }
          else
          {
            fitTxt.push_back("Fit was not performed successfully.");
            fitTxt.push_back(TString::Format("fit status = %d", fitStatus).Data());
          }

          fitTxt.push_back("");
          fitTxt.push_back("Per-bin values (displayed bins only):");
          fitTxt.push_back("binLo  binHi  data_dN_over_NdE  data_err  sim_dN_over_NdE  sim_err  ratio  ratio_err  usedInFit");

          for (size_t i = 0; i < bins.size(); ++i)
          {
            const double ratio  = (i < gyRatio.size()) ? gyRatio[i] : 0.0;
            const double eratio = (i < geyRatio.size()) ? geyRatio[i] : 0.0;
            const bool usedInFit = (ratio > 0.0 && eratio > 0.0 && std::isfinite(ratio) && std::isfinite(eratio));

            fitTxt.push_back(
              TString::Format("%.6g  %.6g  %.10g  %.10g  %.10g  %.10g  %.10g  %.10g  %d",
                              bins[i].lo,
                              bins[i].hi,
                              gyData[i],
                              geyData[i],
                              gySim[i],
                              geySim[i],
                              ratio,
                              eratio,
                              (usedInFit ? 1 : 0)).Data()
            );
          }

          WriteTextFile(
            JoinPath(outYieldDir, "figure30_like_purityCorrectedData_vs_pythiaSignal.txt"),
            fitTxt
          );

          if (gPull) delete gPull;
          delete pTop;
          delete pMid;
          delete pBot;
        };

        auto MakeSamOverlayRawPurityPlot =
          [&]()->void
        {
          if (gApplyPurityCorrectionForUnfolding) return;

          const string samDir = JoinPath(phoDir, "SamOverlay");
          EnsureDir(samDir);

          const string kBoxCutsFile = InputPP();
          const string kBDTFile =
                        "/Users/patsfan753/Desktop/ThesisAnalysis/InputFiles/SamsInputFiles/histsData.root";

          const int nSamBins = 9;
                      const double samPtEdges[nSamBins + 1] = {10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 35.0};

          auto ComputeRawPurityFromCounts =
            [&](double A, double eA,
                double B, double eB,
                double C, double eC,
                double D, double eD,
                double& P, double& eP)->bool
          {
            P = 0.0;
            eP = 0.0;

            if (A <= 0.0 || D <= 0.0) return false;

            double SA = A - B * (C / D);
            if (SA < 0.0) SA = 0.0;
            P = SA / A;

            const double dPdA =  (B * C) / (A * A * D);
            const double dPdB = -(C) / (A * D);
            const double dPdC = -(B) / (A * D);
            const double dPdD =  (B * C) / (A * D * D);

            double var = 0.0;
            var += dPdA * dPdA * eA * eA;
            var += dPdB * dPdB * eB * eB;
            var += dPdC * dPdC * eC * eC;
            var += dPdD * dPdD * eD * eD;

            eP = (var > 0.0) ? std::sqrt(var) : 0.0;
            return true;
          };

          auto IntegrateHistWithOverlap =
            [&](TH1* h, double xLo, double xHi, double& err)->double
          {
            err = 0.0;
            if (!h) return 0.0;

            double sum = 0.0;
            double err2 = 0.0;

            const int nb = h->GetNbinsX();
            for (int ib = 1; ib <= nb; ++ib)
            {
              const double bLo = h->GetXaxis()->GetBinLowEdge(ib);
              const double bHi = h->GetXaxis()->GetBinUpEdge(ib);
              const double ovLo = std::max(xLo, bLo);
              const double ovHi = std::min(xHi, bHi);
              if (ovHi <= ovLo) continue;

              const double width = bHi - bLo;
              if (width <= 0.0) continue;

              const double frac = (ovHi - ovLo) / width;
              const double val = h->GetBinContent(ib);
              double eBin = h->GetBinError(ib);
              if (eBin <= 0.0 && val > 0.0) eBin = std::sqrt(val);

              sum  += frac * val;
              err2 += frac * frac * eBin * eBin;
            }

            err = (err2 > 0.0) ? std::sqrt(err2) : 0.0;
            return sum;
          };

          auto GetObjFromMaybeDir =
            [&](TDirectory* dPrimary, TDirectory* dFallback, const string& name)->TObject*
          {
            TObject* obj = nullptr;
            if (dPrimary) obj = dPrimary->Get(name.c_str());
            if (!obj && dFallback) obj = dFallback->Get(name.c_str());
            return obj;
          };

          TFile* fBox = TFile::Open(kBoxCutsFile.c_str(), "READ");
          if (!fBox || fBox->IsZombie())
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] SamOverlay purity plot skipped: could not open box-cuts file: "
                 << kBoxCutsFile << ANSI_RESET << "\n";
            if (fBox) { fBox->Close(); delete fBox; }
            return;
          }

          TFile* fBDT = TFile::Open(kBDTFile.c_str(), "READ");
          if (!fBDT || fBDT->IsZombie())
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] SamOverlay purity plot skipped: could not open BDT file: "
                 << kBDTFile << ANSI_RESET << "\n";
            if (fBox) { fBox->Close(); delete fBox; }
            if (fBDT) { fBDT->Close(); delete fBDT; }
            return;
          }

          TDirectory* dBox = fBox->GetDirectory(kTriggerPP.c_str());
          if (!dBox) dBox = fBox;

          TDirectory* dBDTTop = fBDT;
          TDirectory* dBDTTrig = fBDT->GetDirectory(kTriggerPP.c_str());

          std::vector<double> xBox, exBox, yBox, eyBox;
          std::vector<double> xBDT, exBDT, yBDT, eyBDT;

          xBox.reserve(nSamBins);
          exBox.reserve(nSamBins);
          yBox.reserve(nSamBins);
          eyBox.reserve(nSamBins);

          xBDT.reserve(nSamBins);
          exBDT.reserve(nSamBins);
          yBDT.reserve(nSamBins);
          eyBDT.reserve(nSamBins);

          TH1* hBDTA = dynamic_cast<TH1*>(GetObjFromMaybeDir(dBDTTop, dBDTTrig, "hclusterptabcd0"));
          TH1* hBDTB = dynamic_cast<TH1*>(GetObjFromMaybeDir(dBDTTop, dBDTTrig, "hclusterptabcd1"));
          TH1* hBDTC = dynamic_cast<TH1*>(GetObjFromMaybeDir(dBDTTop, dBDTTrig, "hclusterptabcd2"));
          TH1* hBDTD = dynamic_cast<TH1*>(GetObjFromMaybeDir(dBDTTop, dBDTTrig, "hclusterptabcd3"));

          if (!hBDTA || !hBDTB || !hBDTC || !hBDTD)
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] SamOverlay purity plot skipped: missing one or more BDT histograms hclusterptabcd{0,1,2,3} in "
                 << kBDTFile << ANSI_RESET << "\n";
            fBox->Close();
            fBDT->Close();
            delete fBox;
            delete fBDT;
            return;
          }

          for (int i = 0; i < nSamBins; ++i)
          {
            const int ptLo = (int)std::llround(samPtEdges[i]);
            const int ptHi = (int)std::llround(samPtEdges[i + 1]);
            const double ptCtr = 0.5 * (samPtEdges[i] + samPtEdges[i + 1]);
            const double ptErr = 0.5 * (samPtEdges[i + 1] - samPtEdges[i]);

            {
              const string hAName = TString::Format("h_isIsolated_isTight_pT_%d_%d", ptLo, ptHi).Data();
              const string hBName = TString::Format("h_notIsolated_isTight_pT_%d_%d", ptLo, ptHi).Data();
              const string hCName = TString::Format("h_isIsolated_notTight_pT_%d_%d", ptLo, ptHi).Data();
              const string hDName = TString::Format("h_notIsolated_notTight_pT_%d_%d", ptLo, ptHi).Data();

              TH1* hA = dynamic_cast<TH1*>(dBox->Get(hAName.c_str()));
              TH1* hB = dynamic_cast<TH1*>(dBox->Get(hBName.c_str()));
              TH1* hC = dynamic_cast<TH1*>(dBox->Get(hCName.c_str()));
              TH1* hD = dynamic_cast<TH1*>(dBox->Get(hDName.c_str()));

              if (hA && hB && hC && hD)
              {
                const double A = hA->GetBinContent(1);
                const double B = hB->GetBinContent(1);
                const double C = hC->GetBinContent(1);
                const double D = hD->GetBinContent(1);

                double eA = hA->GetBinError(1);
                double eB = hB->GetBinError(1);
                double eC = hC->GetBinError(1);
                double eD = hD->GetBinError(1);

                if (eA <= 0.0 && A > 0.0) eA = std::sqrt(A);
                if (eB <= 0.0 && B > 0.0) eB = std::sqrt(B);
                if (eC <= 0.0 && C > 0.0) eC = std::sqrt(C);
                if (eD <= 0.0 && D > 0.0) eD = std::sqrt(D);

                double P = 0.0;
                double eP = 0.0;
                if (ComputeRawPurityFromCounts(A, eA, B, eB, C, eC, D, eD, P, eP))
                {
                  xBox.push_back(ptCtr);
                  exBox.push_back(ptErr);
                  yBox.push_back(P);
                  eyBox.push_back(eP);
                }
              }
              else
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] SamOverlay box-cuts purity: missing one or more ABCD histograms for pT bin "
                     << ptLo << "-" << ptHi << " GeV in " << kBoxCutsFile
                     << ANSI_RESET << "\n";
              }
            }

            {
              double eA = 0.0, eB = 0.0, eC = 0.0, eD = 0.0;
              const double A = IntegrateHistWithOverlap(hBDTA, samPtEdges[i], samPtEdges[i + 1], eA);
              const double B = IntegrateHistWithOverlap(hBDTB, samPtEdges[i], samPtEdges[i + 1], eB);
              const double C = IntegrateHistWithOverlap(hBDTC, samPtEdges[i], samPtEdges[i + 1], eC);
              const double D = IntegrateHistWithOverlap(hBDTD, samPtEdges[i], samPtEdges[i + 1], eD);

              double P = 0.0;
              double eP = 0.0;
              if (ComputeRawPurityFromCounts(A, eA, B, eB, C, eC, D, eD, P, eP))
              {
                xBDT.push_back(ptCtr);
                exBDT.push_back(ptErr);
                yBDT.push_back(P);
                eyBDT.push_back(eP);
              }
              else
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] SamOverlay BDT purity: invalid raw ABCD counts for pT bin "
                     << ptLo << "-" << ptHi << " GeV in " << kBDTFile
                     << ANSI_RESET << "\n";
              }
            }
          }

          if (xBox.empty() || xBDT.empty())
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] SamOverlay purity plot skipped: no valid overlay points were built."
                 << ANSI_RESET << "\n";
            fBox->Close();
            fBDT->Close();
            delete fBox;
            delete fBDT;
            return;
          }

          TGraphErrors gBDT((int)xBDT.size(), &xBDT[0], &yBDT[0], &exBDT[0], &eyBDT[0]);
          TGraphErrors gBox((int)xBox.size(), &xBox[0], &yBox[0], &exBox[0], &eyBox[0]);

          gBDT.SetTitle("");
          gBDT.SetMarkerStyle(20);
          gBDT.SetMarkerSize(1.20);
          gBDT.SetMarkerColor(kBlue + 1);
          gBDT.SetLineColor(kBlue + 1);
          gBDT.SetLineWidth(2);

          gBox.SetTitle("");
          gBox.SetMarkerStyle(20);
          gBox.SetMarkerSize(1.20);
          gBox.SetMarkerColor(kRed + 1);
          gBox.SetLineColor(kRed + 1);
          gBox.SetLineWidth(2);

          TCanvas c("c_samOverlay_purityRaw", "c_samOverlay_purityRaw", 900, 700);
          ApplyCanvasMargins1D(c);

          TH1F hFrame("hSamOverlayPurityFrame", "", 100, samPtEdges[0], samPtEdges[nSamBins]);
          hFrame.SetDirectory(nullptr);
          hFrame.SetStats(0);
          hFrame.SetTitle("");
          hFrame.SetMinimum(0.0);
          hFrame.SetMaximum(1.05);
          hFrame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
          hFrame.GetYaxis()->SetTitle("Purity (raw ABCD)");
          hFrame.Draw();

          gBox.Draw("PE same");
          gBDT.Draw("PE same");

          TLegend leg(0.58, 0.76, 0.88, 0.88);
          leg.SetBorderSize(0);
          leg.SetFillStyle(0);
          leg.SetTextFont(42);
          leg.SetTextSize(0.035);
          leg.AddEntry(&gBDT, "BDT tight-ID", "pe");
          leg.AddEntry(&gBox, "Box-cuts tight-ID", "pe");
          leg.Draw();

          TLatex tTitle;
          tTitle.SetNDC(true);
          tTitle.SetTextFont(42);
          tTitle.SetTextAlign(23);
          tTitle.SetTextSize(0.045);
          tTitle.DrawLatex(0.50, 0.96, "Purity Overlay, BDT vs Box Cuts");

          TLatex tInfo;
          tInfo.SetNDC(true);
          tInfo.SetTextFont(42);
          tInfo.SetTextAlign(33);
          tInfo.SetTextSize(0.038);
          tInfo.DrawLatex(0.92, 0.45, "Iso cone: 0.03");
          tInfo.DrawLatex(0.92, 0.35, "|v_{z}| < 30 cm");

          SaveCanvas(c, JoinPath(samDir, "purity_raw_overlay_BDT_vs_BoxCuts.png"));

          fBox->Close();
          fBDT->Close();
          delete fBox;
          delete fBDT;
        };

        MakeSamOverlayRawPurityPlot();

        // -----------------------------------------------------------------
        // 2x2 ABCD counts comparison panel: last 2 pT bins, Box vs BDT
        // rows = pT bin (24-26, 26-35), cols = Box cuts | BDT (Sam)
        // -----------------------------------------------------------------
        auto MakeABCDComparisonPanel = [&]()->void
        {
          const string kBDTFile2 =
                        "/Users/patsfan753/Desktop/ThesisAnalysis/InputFiles/SamsInputFiles/histsData.root";

          TFile* fBDT2 = TFile::Open(kBDTFile2.c_str(), "READ");
          if (!fBDT2 || fBDT2->IsZombie())
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] ABCDComparePanel skipped: could not open BDT file: "
                 << kBDTFile2 << ANSI_RESET << "\n";
            if (fBDT2) { fBDT2->Close(); delete fBDT2; }
            return;
          }

          TDirectory* dBDTTop2  = fBDT2;
          TDirectory* dBDTTrig2 = fBDT2->GetDirectory(kTriggerPP.c_str());

          TH1* hBDTA2 = dynamic_cast<TH1*>(dBDTTop2->Get("hclusterptabcd0"));
          if (!hBDTA2 && dBDTTrig2) hBDTA2 = dynamic_cast<TH1*>(dBDTTrig2->Get("hclusterptabcd0"));
          TH1* hBDTB2 = dynamic_cast<TH1*>(dBDTTop2->Get("hclusterptabcd1"));
          if (!hBDTB2 && dBDTTrig2) hBDTB2 = dynamic_cast<TH1*>(dBDTTrig2->Get("hclusterptabcd1"));
          TH1* hBDTC2 = dynamic_cast<TH1*>(dBDTTop2->Get("hclusterptabcd2"));
          if (!hBDTC2 && dBDTTrig2) hBDTC2 = dynamic_cast<TH1*>(dBDTTrig2->Get("hclusterptabcd2"));
          TH1* hBDTD2 = dynamic_cast<TH1*>(dBDTTop2->Get("hclusterptabcd3"));
          if (!hBDTD2 && dBDTTrig2) hBDTD2 = dynamic_cast<TH1*>(dBDTTrig2->Get("hclusterptabcd3"));

          if (!hBDTA2 || !hBDTB2 || !hBDTC2 || !hBDTD2)
          {
            cout << ANSI_BOLD_YEL
                 << "[WARN] ABCDComparePanel skipped: missing hclusterptabcd{0-3} in "
                 << kBDTFile2 << ANSI_RESET << "\n";
            fBDT2->Close(); delete fBDT2;
            return;
          }

          auto IntegrateSlice = [](TH1* h, double xLo, double xHi)->double
          {
            if (!h) return 0.0;
            double sum = 0.0;
            const int nb = h->GetNbinsX();
            for (int ib = 1; ib <= nb; ++ib)
            {
              const double bLo = h->GetXaxis()->GetBinLowEdge(ib);
              const double bHi = h->GetXaxis()->GetBinUpEdge(ib);
              const double ovLo = std::max(xLo, bLo);
              const double ovHi = std::min(xHi, bHi);
              if (ovHi <= ovLo) continue;
              const double width = bHi - bLo;
              if (width <= 0.0) continue;
              sum += (ovHi - ovLo) / width * h->GetBinContent(ib);
            }
            return sum;
          };

          // pT bins to display: last two from kPtEdges = {...,24,26,35}
          struct PanelBin { int lo; int hi; const char* suf; };
          const PanelBin panelBins[2] = {
            {24, 26, "_pT_24_26"},
            {26, 35, "_pT_26_35"}
          };

          const char* xLabelsABCD[4]  = {"N_{A}", "N_{B}", "N_{C}", "N_{D}"};
          const int   colorsABCD[4]   = {kGreen+2, kRed+1, kAzure+1, kMagenta+1};

          // Canvas: 2 columns x 2 rows  (col1=Box, col2=BDT; row1=24-26, row2=26-35)
          TCanvas cPanel("c_abcd_compare_panel","c_abcd_compare_panel", 1800, 1600);
          cPanel.Divide(2, 2, 0.002, 0.002);

          vector<TObject*> keepAlive2;
          keepAlive2.reserve(2 * 2 * 5);

          auto DrawABCDBarPad = [&](int padIdx,
                                    double A, double B, double Cc, double D,
                                    const char* colTitle,
                                    int ptLo, int ptHi)
          {
            cPanel.cd(padIdx);
            if (!gPad) return;

            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.14);
            gPad->SetBottomMargin(0.22);
            gPad->SetTicks(1,1);

            const double vals[4] = {A, B, Cc, D};
            double ymax = 0.0;
            for (int ib = 0; ib < 4; ++ib) ymax = std::max(ymax, vals[ib]);
            const double yMaxPlot = (ymax > 0.0) ? (1.35 * ymax) : 1.0;

            TH1F* hAxis = new TH1F(
              TString::Format("h_cmpAxis_p%d", padIdx).Data(),
              "", 4, 0.5, 4.5
            );
            hAxis->SetDirectory(nullptr);
            hAxis->SetStats(0);
            hAxis->SetMinimum(0.0);
            hAxis->SetMaximum(yMaxPlot);
            for (int ib = 1; ib <= 4; ++ib) hAxis->GetXaxis()->SetBinLabel(ib, xLabelsABCD[ib-1]);
            hAxis->GetYaxis()->SetTitle("Counts");
            hAxis->GetXaxis()->SetTitle("");
            hAxis->GetXaxis()->LabelsOption("h");
            hAxis->GetXaxis()->SetLabelSize(0.070);
            hAxis->GetXaxis()->SetLabelOffset(0.010);
            hAxis->GetYaxis()->SetTitleSize(0.054);
            hAxis->GetYaxis()->SetTitleOffset(1.28);
            hAxis->GetYaxis()->SetLabelSize(0.046);
            hAxis->SetLineColor(1);
            hAxis->SetLineWidth(2);
            hAxis->SetFillStyle(0);
            hAxis->Draw("hist");
            keepAlive2.push_back(hAxis);

            for (int ib = 1; ib <= 4; ++ib)
            {
              TH1F* hb = new TH1F(
                TString::Format("h_cmpBar_p%d_b%d", padIdx, ib).Data(),
                "", 4, 0.5, 4.5
              );
              hb->SetDirectory(nullptr);
              hb->SetStats(0);
              hb->SetBinContent(ib, vals[ib-1]);
              hb->SetFillStyle(1001);
              hb->SetFillColor(colorsABCD[ib-1]);
              hb->SetLineColor(1);
              hb->SetLineWidth(2);
              hb->SetBarWidth(0.90);
              hb->SetBarOffset(0.05);
              hb->Draw("BAR SAME");
              keepAlive2.push_back(hb);
            }

            TLatex t;
            t.SetTextFont(42);
            t.SetTextAlign(22);
            t.SetTextSize(0.052);
            for (int ib = 1; ib <= 4; ++ib)
            {
              const double y = vals[ib-1];
              if (y <= 0.0) continue;
              const double x = hAxis->GetXaxis()->GetBinCenter(ib);
              const double yText = std::min(y + 0.025*yMaxPlot, 0.93*yMaxPlot);
              t.DrawLatex(x, yText, TString::Format("%.0f", y).Data());
            }

            TLatex tt;
            tt.SetTextFont(42);
            tt.SetNDC();
            tt.SetTextAlign(23);
            tt.SetTextSize(0.058);
            tt.DrawLatex(0.50, 0.965,
              TString::Format("p_{T}^{#gamma}: %d-%d GeV  [%s]", ptLo, ptHi, colTitle).Data()
            );
          };

          // Row 1: pT 24-26  |  Row 2: pT 26-35
          // Pad ordering: Divide(2,2) -> pad 1=top-left, 2=top-right, 3=bot-left, 4=bot-right
          for (int row = 0; row < 2; ++row)
          {
            const PanelBin& pb = panelBins[row];
            const string sufStr = pb.suf;

            // Box cuts column (pad 2*row+1)
            {
              const double A  = Read1BinCount(dsData, "h_isIsolated_isTight"  + sufStr);
              const double B  = Read1BinCount(dsData, "h_notIsolated_isTight" + sufStr);
              const double Cc = Read1BinCount(dsData, "h_isIsolated_notTight" + sufStr);
              const double D  = Read1BinCount(dsData, "h_notIsolated_notTight"+ sufStr);
              DrawABCDBarPad(2*row + 1, A, B, Cc, D, "Box cuts", pb.lo, pb.hi);
            }

            // BDT column (pad 2*row+2)
            {
              const double A  = IntegrateSlice(hBDTA2, (double)pb.lo, (double)pb.hi);
              const double B  = IntegrateSlice(hBDTB2, (double)pb.lo, (double)pb.hi);
              const double Cc = IntegrateSlice(hBDTC2, (double)pb.lo, (double)pb.hi);
              const double D  = IntegrateSlice(hBDTD2, (double)pb.lo, (double)pb.hi);
              DrawABCDBarPad(2*row + 2, A, B, Cc, D, "BDT", pb.lo, pb.hi);
            }
          }

          cPanel.cd(0);

          SaveCanvas(cPanel, JoinPath(JoinPath(phoDir, "SamOverlay"), "abcd_counts_compare_BoxCuts_vs_BDT.png"));

          for (auto* obj : keepAlive2) delete obj;
          fBDT2->Close();
          delete fBDT2;
        };

        MakeABCDComparisonPanel();

        const string phoRecoName   = "h_unfoldRecoPho_pTgamma";
        const string phoTruthName  = "h_unfoldTruthPho_pTgamma";
        const string phoRespName   = "h2_unfoldResponsePho_pTgamma";
        const string phoFakesName  = "h_unfoldRecoPhoFakes_pTgamma";       // SIM-only bookkeeping
        const string phoMissesName = "h_unfoldTruthPhoMisses_pTgamma";     // SIM-only bookkeeping

        TH1* hPhoRecoData_in        = GetObj<TH1>(dsData, phoRecoName,   true, true, true);
        TH1* hPhoRecoSim_in         = GetObj<TH1>(dsSim,  phoRecoName,   true, true, true);
        TH1* hPhoTruthSim_in        = GetObj<TH1>(dsSim,  phoTruthName,  true, true, true);
        TH2* hPhoRespSim_in         = GetObj<TH2>(dsSim,  phoRespName,   true, true, true);
        TH1* hPhoRecoFakesSim_in    = GetObj<TH1>(dsSim,  phoFakesName,  true, true, true);
        TH1* hPhoTruthMissesSim_in  = GetObj<TH1>(dsSim,  phoMissesName, true, true, true);

        auto MissingWhy = [&](Dataset& ds, const string& relName)->string
        {
          const string fp = FullPath(ds, relName);
          auto it = ds.missingReason.find(fp);
          if (it == ds.missingReason.end()) return "";
          return it->second;
        };

        auto PrintGet = [&](Dataset& ds, const string& relName, TObject* obj)
        {
          const string fp = FullPath(ds, relName);
          cout << "  [GET] " << ds.label << " :: " << fp << " : ";

          if (!obj)
          {
            const string why = MissingWhy(ds, relName);
            cout << ANSI_BOLD_RED << "MISSING" << ANSI_RESET;
            if (!why.empty()) cout << " (" << why << ")";
            cout << "\n";
            return;
          }

          cout << ANSI_BOLD_GRN << "FOUND" << ANSI_RESET << " (" << obj->ClassName() << ")";
          if (auto* h2 = dynamic_cast<TH2*>(obj))
          {
            cout << "  entries=" << h2->GetEntries()
                 << "  nbinsX=" << h2->GetNbinsX()
                 << "  nbinsY=" << h2->GetNbinsY();
          }
          else if (auto* h1 = dynamic_cast<TH1*>(obj))
          {
            cout << "  entries=" << h1->GetEntries()
                 << "  nbinsX=" << h1->GetNbinsX();
          }
          cout << "\n";
        };

        auto DumpTH1Summary = [&](const string& tag, TH1* h)->void
        {
          cout << ANSI_BOLD_CYN << "[DEBUG TH1] " << tag << ANSI_RESET << "\n";
          if (!h)
          {
            cout << "  <null>\n";
            return;
          }

          cout << "  name=" << h->GetName()
               << "  title=" << h->GetTitle()
               << "  nbins=" << h->GetNbinsX()
               << "  entries=" << h->GetEntries()
               << "  integral(width)=" << h->Integral(1, h->GetNbinsX(), "width")
               << "  integral=" << h->Integral(1, h->GetNbinsX())
               << "\n";

          for (int ib = 0; ib <= h->GetNbinsX() + 1; ++ib)
          {
            const double lo  = (ib >= 1 && ib <= h->GetNbinsX()) ? h->GetXaxis()->GetBinLowEdge(ib) : -999.0;
            const double hi  = (ib >= 1 && ib <= h->GetNbinsX()) ? h->GetXaxis()->GetBinUpEdge(ib)  : -999.0;
            const double val = h->GetBinContent(ib);
            const double err = h->GetBinError(ib);

            if (ib == 0)
            {
              cout << "    UF"
                   << "  content=" << val
                   << "  error="   << err << "\n";
            }
            else if (ib == h->GetNbinsX() + 1)
            {
              cout << "    OF"
                   << "  content=" << val
                   << "  error="   << err << "\n";
            }
            else
            {
              cout << "    bin " << ib
                   << "  [" << lo << "," << hi << "]"
                   << "  content=" << val
                   << "  error="   << err << "\n";
            }
          }
        };

        auto DumpTH2Summary = [&](const string& tag, TH2* h)->void
        {
          cout << ANSI_BOLD_CYN << "[DEBUG TH2] " << tag << ANSI_RESET << "\n";
          if (!h)
          {
            cout << "  <null>\n";
            return;
          }

          cout << "  name=" << h->GetName()
               << "  title=" << h->GetTitle()
               << "  nbinsX=" << h->GetNbinsX()
               << "  nbinsY=" << h->GetNbinsY()
               << "  entries=" << h->GetEntries()
               << "  integral=" << h->Integral(1, h->GetNbinsX(), 1, h->GetNbinsY())
               << "\n";

          for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
          {
            double rowInt = 0.0;
            double rowIntWidth = 0.0;
            for (int iy = 1; iy <= h->GetNbinsY(); ++iy)
            {
              rowInt += h->GetBinContent(ix, iy);
              rowIntWidth += h->GetBinContent(ix, iy) * h->GetYaxis()->GetBinWidth(iy);
            }

            cout << "    x-bin " << ix
                 << "  [" << h->GetXaxis()->GetBinLowEdge(ix)
                 << ","   << h->GetXaxis()->GetBinUpEdge(ix) << "]"
                 << "  rowIntegral=" << rowInt
                 << "  rowIntegral(widthY)=" << rowIntWidth
                 << "\n";
          }
        };

        auto DumpProjectionYSummary = [&](const string& tag, TH2* h2, int ix)->void
        {
          cout << ANSI_BOLD_CYN << "[DEBUG PROJY] " << tag << ANSI_RESET << "\n";
          if (!h2)
          {
            cout << "  source TH2 is null\n";
            return;
          }
          if (ix < 1 || ix > h2->GetNbinsX())
          {
            cout << "  requested ix=" << ix << " is out of range 1.." << h2->GetNbinsX() << "\n";
            return;
          }

          TH1D* hp = h2->ProjectionY(
            TString::Format("hDebugProjY_%s_ix%d", h2->GetName(), ix).Data(),
            ix, ix, "e"
          );
          if (!hp)
          {
            cout << "  ProjectionY returned null\n";
            return;
          }

          hp->SetDirectory(nullptr);
          EnsureSumw2(hp);

          cout << "  x-bin " << ix
               << "  [" << h2->GetXaxis()->GetBinLowEdge(ix)
               << ","   << h2->GetXaxis()->GetBinUpEdge(ix) << "]"
               << "  projIntegral=" << hp->Integral(1, hp->GetNbinsX())
               << "  projIntegral(width)=" << hp->Integral(1, hp->GetNbinsX(), "width")
               << "\n";

          for (int iy = 0; iy <= hp->GetNbinsX() + 1; ++iy)
          {
            if (iy == 0)
            {
              cout << "    UF"
                   << "  content=" << hp->GetBinContent(iy)
                   << "  error="   << hp->GetBinError(iy) << "\n";
            }
            else if (iy == hp->GetNbinsX() + 1)
            {
              cout << "    OF"
                   << "  content=" << hp->GetBinContent(iy)
                   << "  error="   << hp->GetBinError(iy) << "\n";
            }
            else
            {
              cout << "    y-bin " << iy
                   << "  [" << hp->GetXaxis()->GetBinLowEdge(iy)
                   << ","   << hp->GetXaxis()->GetBinUpEdge(iy) << "]"
                   << "  content=" << hp->GetBinContent(iy)
                   << "  error="   << hp->GetBinError(iy) << "\n";
            }
          }

          delete hp;
        };

        struct IterScanSummary
        {
            vector<double> xIt;
            vector<double> exIt;
            vector<double> yRelStat;
            vector<double> eyRelStat;
            vector<double> yRelChange;
            vector<double> eyRelChange;
            vector<double> yQuad;
            vector<double> eyQuad;
            int bestIt = -1;
            double bestQuad = std::numeric_limits<double>::max();
        };

        auto BuildPhotonIterScan = [&](RooUnfoldResponse& resp,
                                         TH1* hRecoMeasured,
                                         TH1* hTruthTemplate,
                                         int kMaxIt,
                                         int nToys)->IterScanSummary
        {
            IterScanSummary s;
            if (!hRecoMeasured || !hTruthTemplate) return s;

            vector<int> selTruthBins;
            {
              const auto& analysisBins = UnfoldAnalysisRecoPtBins();
              TAxis* axTruth = hTruthTemplate->GetXaxis();

              for (const auto& b : analysisBins)
              {
                int ibMatch = -1;
                for (int ib = 1; ib <= axTruth->GetNbins(); ++ib)
                {
                  const double lo = axTruth->GetBinLowEdge(ib);
                  const double hi = axTruth->GetBinUpEdge(ib);
                  if (std::fabs(lo - (double)b.lo) < 1e-9 &&
                      std::fabs(hi - (double)b.hi) < 1e-9)
                  {
                    ibMatch = ib;
                    break;
                  }
                }

                if (ibMatch > 0) selTruthBins.push_back(ibMatch);
              }
            }

            if (selTruthBins.empty()) return s;

            TH1* hPrev = CloneTH1(hTruthTemplate, "hPhoIter0Baseline");
            if (!hPrev) return s;

            hPrev->SetDirectory(nullptr);
            EnsureSumw2(hPrev);
            hPrev->Reset("ICES");

            for (int ib = 1; ib <= hPrev->GetNbinsX(); ++ib)
            {
              const double x = hPrev->GetXaxis()->GetBinCenter(ib);
              const int iReco = hRecoMeasured->GetXaxis()->FindBin(x);
              if (iReco < 1 || iReco > hRecoMeasured->GetNbinsX()) continue;

              hPrev->SetBinContent(ib, hRecoMeasured->GetBinContent(iReco));
              hPrev->SetBinError  (ib, hRecoMeasured->GetBinError  (iReco));
            }

            for (int it = 1; it <= kMaxIt; ++it)
            {
              RooUnfoldBayes uIt(&resp, hRecoMeasured, it);
              uIt.SetVerbose(0);
              uIt.SetNToys(nToys);

              TH1* hIt = nullptr;
              if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
              hIt = uIt.Hreco(RooUnfold::kCovToy);
              const auto cov = uIt.Ereco(RooUnfold::kCovToy);
              if (gSystem) gSystem->RedirectOutput(0);
              if (!hIt) continue;

              hIt->SetDirectory(nullptr);
              EnsureSumw2(hIt);

              double sumV = 0.0;
              for (const int ib : selTruthBins)
              {
                if (ib < 1 || ib > hIt->GetNbinsX()) continue;
                sumV += hIt->GetBinContent(ib);
              }

              double sumCov = 0.0;
              const int nCovRows = cov.GetNrows();
              const int nCovCols = cov.GetNcols();

              for (size_t i = 0; i < selTruthBins.size(); ++i)
              {
                const int ii = selTruthBins[i] - 1;
                if (ii < 0 || ii >= nCovRows) continue;

                for (size_t j = 0; j < selTruthBins.size(); ++j)
                {
                  const int jj = selTruthBins[j] - 1;
                  if (jj < 0 || jj >= nCovCols) continue;
                  sumCov += cov(ii, jj);
                }
              }

              double relStat = 0.0;
              if (sumV > 0.0) relStat = std::sqrt(std::max(0.0, sumCov)) / sumV;

              double relDev = 0.0;
              {
                double num2 = 0.0;
                double den2 = 0.0;

                for (const int ib : selTruthBins)
                {
                  if (ib < 1 || ib > hIt->GetNbinsX()) continue;

                  const double v  = hIt->GetBinContent(ib);
                  const double vp = hPrev->GetBinContent(ib);
                  const double d  = v - vp;

                  num2 += d * d;
                  den2 += v * v;
                }

                if (den2 > 0.0) relDev = std::sqrt(num2 / den2);
              }

              const double quad = std::sqrt(relStat * relStat + relDev * relDev);

              s.xIt.push_back((double)it);
              s.exIt.push_back(0.0);
              s.yRelStat.push_back(relStat);
              s.eyRelStat.push_back(0.0);
              s.yRelChange.push_back(relDev);
              s.eyRelChange.push_back(0.0);
              s.yQuad.push_back(quad);
              s.eyQuad.push_back(0.0);

              if (quad < s.bestQuad)
              {
                s.bestQuad = quad;
                s.bestIt = it;
              }

              delete hPrev;
              hPrev = hIt;
            }

            if (hPrev) delete hPrev;

            return s;
          };

          auto BuildXJIterScan = [&](RooUnfoldResponse& resp,
                                       TH1* hMeasuredGlob,
                                       TH1* hTruthTemplateGlob,
                                       TH2* hMeasured2D,
                                       TH2* hTruthTemplate2D,
                                       int kMaxIt,
                                       int nToys,
                                       const string& nameSuffix)->IterScanSummary
          {
              IterScanSummary s;
              if (!hMeasuredGlob || !hTruthTemplateGlob || !hMeasured2D || !hTruthTemplate2D) return s;

              const int nGlobMeasured2D = (hMeasured2D->GetNbinsX() + 2) * (hMeasured2D->GetNbinsY() + 2);
              const int nGlobTruth2D    = (hTruthTemplate2D->GetNbinsX() + 2) * (hTruthTemplate2D->GetNbinsY() + 2);

              if (hMeasuredGlob->GetNbinsX() != nGlobMeasured2D) return s;
              if (hTruthTemplateGlob->GetNbinsX() != nGlobTruth2D) return s;

              vector<int> selTruthXBins;
              vector<int> selTruthGlobalBins;

              {
                const auto& analysisBins = UnfoldAnalysisRecoPtBins();
                TAxis* axTruth = hTruthTemplate2D->GetXaxis();

                for (const auto& b : analysisBins)
                {
                  int ixMatch = -1;
                  for (int ix = 1; ix <= axTruth->GetNbins(); ++ix)
                  {
                    const double lo = axTruth->GetBinLowEdge(ix);
                    const double hi = axTruth->GetBinUpEdge(ix);
                    if (std::fabs(lo - (double)b.lo) < 1e-9 &&
                        std::fabs(hi - (double)b.hi) < 1e-9)
                    {
                      ixMatch = ix;
                      break;
                    }
                  }

                  if (ixMatch > 0) selTruthXBins.push_back(ixMatch);
                }

                const int nyTruth = hTruthTemplate2D->GetNbinsY();
                for (const int ix : selTruthXBins)
                {
                  for (int iy = 1; iy <= nyTruth; ++iy)
                  {
                    selTruthGlobalBins.push_back(hTruthTemplate2D->GetBin(ix, iy));
                  }
                }
              }

              if (selTruthXBins.empty() || selTruthGlobalBins.empty()) return s;

              TH1* hPrev = CloneTH1(hTruthTemplateGlob, "hXJIter0Baseline");
              if (!hPrev) return s;

              hPrev->SetDirectory(nullptr);
              EnsureSumw2(hPrev);
              hPrev->Reset("ICES");

              const int nbPrevFill = std::min(hPrev->GetNbinsX(), hMeasuredGlob->GetNbinsX());
              for (int ib = 1; ib <= nbPrevFill; ++ib)
              {
                hPrev->SetBinContent(ib, hMeasuredGlob->GetBinContent(ib));
                hPrev->SetBinError  (ib, hMeasuredGlob->GetBinError  (ib));
              }

              for (int it = 1; it <= kMaxIt; ++it)
              {
                RooUnfoldBayes u(&resp, hMeasuredGlob, it);
                u.SetVerbose(0);
                u.SetNToys(nToys);

                TH1* hCurr = nullptr;
                if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
                hCurr = u.Hreco(RooUnfold::kCovToy);
                const auto cov = u.Ereco(RooUnfold::kCovToy);
                if (gSystem) gSystem->RedirectOutput(0);
                if (!hCurr) continue;

                hCurr->SetDirectory(nullptr);
                EnsureSumw2(hCurr);

                if (hCurr->GetNbinsX() != hTruthTemplateGlob->GetNbinsX())
                {
                  delete hCurr;
                  continue;
                }

                TH2* hCurr2D = UnflattenGlobalToTH2(
                  hCurr,
                  hTruthTemplate2D,
                  TString::Format("hXJIterCurr2D_%s_it%d", nameSuffix.c_str(), it).Data()
                );
                TH2* hPrev2D = UnflattenGlobalToTH2(
                  hPrev,
                  hTruthTemplate2D,
                  TString::Format("hXJIterPrev2D_%s_it%d", nameSuffix.c_str(), it).Data()
                );

                if (!hCurr2D || !hPrev2D)
                {
                  if (hCurr2D) delete hCurr2D;
                  if (hPrev2D) delete hPrev2D;
                  delete hPrev;
                  hPrev = hCurr;
                  continue;
                }

                hCurr2D->SetDirectory(nullptr);
                hPrev2D->SetDirectory(nullptr);

                double sumV  = 0.0;
                double sumV2 = 0.0;
                double sumD2 = 0.0;

                const int nyTruth = hCurr2D->GetNbinsY();

                for (const int ix : selTruthXBins)
                {
                  if (ix < 1 || ix > hCurr2D->GetNbinsX()) continue;

                  for (int iy = 1; iy <= nyTruth; ++iy)
                  {
                    const double v  = hCurr2D->GetBinContent(ix, iy);
                    const double vp = hPrev2D->GetBinContent(ix, iy);
                    const double d  = v - vp;

                    if (v == 0.0 && vp == 0.0) continue;

                    sumV  += v;
                    sumV2 += v * v;
                    sumD2 += d * d;
                  }
                }

                double sumCov = 0.0;
                const int nCovRows = cov.GetNrows();
                const int nCovCols = cov.GetNcols();

                for (size_t i = 0; i < selTruthGlobalBins.size(); ++i)
                {
                  const int ii = selTruthGlobalBins[i];
                  if (ii < 0 || ii >= nCovRows) continue;

                  for (size_t j = 0; j < selTruthGlobalBins.size(); ++j)
                  {
                    const int jj = selTruthGlobalBins[j];
                    if (jj < 0 || jj >= nCovCols) continue;
                    sumCov += cov(ii, jj);
                  }
                }

                const double relStat = (sumV > 0.0) ? (std::sqrt(std::max(0.0, sumCov)) / sumV) : 0.0;
                const double relChg  = (sumV2 > 0.0) ? std::sqrt(sumD2 / sumV2) : 0.0;

                s.xIt.push_back((double)it);
                s.exIt.push_back(0.0);
                s.yRelStat.push_back(relStat);
                s.eyRelStat.push_back(0.0);
                s.yRelChange.push_back(relChg);
                s.eyRelChange.push_back(0.0);
                s.yQuad.push_back(std::sqrt(relStat * relStat + relChg * relChg));
                s.eyQuad.push_back(0.0);

                if (s.yQuad.back() < s.bestQuad)
                {
                  s.bestQuad = s.yQuad.back();
                  s.bestIt = it;
                }

                delete hCurr2D;
                delete hPrev2D;
                delete hPrev;
                hPrev = hCurr;
              }

              if (hPrev) delete hPrev;

              return s;
        };
     
        cout << "\n  [5I] Photon unfolding inputs:\n";
        PrintGet(dsData, phoRecoName,   hPhoRecoData_in);
        PrintGet(dsSim,  phoRecoName,   hPhoRecoSim_in);
        PrintGet(dsSim,  phoTruthName,  hPhoTruthSim_in);
        PrintGet(dsSim,  phoRespName,   hPhoRespSim_in);
        PrintGet(dsSim,  phoFakesName,  hPhoRecoFakesSim_in);
        PrintGet(dsSim,  phoMissesName, hPhoTruthMissesSim_in);

        if (!hPhoRecoData_in || !hPhoRecoSim_in || !hPhoTruthSim_in || !hPhoRespSim_in)
        {
          cout << ANSI_BOLD_RED
               << "[ERROR] Missing one or more photon unfolding inputs.\n"
               << "        Need (DATA) h_unfoldRecoPho_pTgamma and (SIM) h_unfoldRecoPho_pTgamma, h_unfoldTruthPho_pTgamma, h2_unfoldResponsePho_pTgamma.\n"
               << "        Aborting RooUnfold pipeline."
               << ANSI_RESET << "\n";
          return;
        }

        TH1* hPhoRecoData  = CloneTH1(hPhoRecoData_in,  "hPhoRecoData");
        TH1* hPhoRecoSim   = CloneTH1(hPhoRecoSim_in,   "hPhoRecoSim");
        TH1* hPhoTruthSim  = CloneTH1(hPhoTruthSim_in,  "hPhoTruthSim");
        TH2* hPhoRespSim   = CloneTH2(hPhoRespSim_in,   "hPhoRespSim_truthVsReco");

        EnsureSumw2(hPhoRecoData);
        EnsureSumw2(hPhoRecoSim);
        EnsureSumw2(hPhoTruthSim);
        EnsureSumw2(hPhoRespSim);

        if (gApplyPurityCorrectionForUnfolding)
        {
          cout << "  [BUILD] DATA :: purity-corrected photon reco input from ABCD counts\n";
          if (!ApplyPurityCorrectionToRecoPhotonHist(hPhoRecoData))
          {
            cout << ANSI_BOLD_RED
                 << "[ERROR] Failed to build purity-corrected photon reco input. Aborting RooUnfold pipeline."
                 << ANSI_RESET << "\n";
            delete hPhoRecoData;
            delete hPhoRecoSim;
            delete hPhoTruthSim;
            delete hPhoRespSim;
            return;
          }
      }

      if (gApplyPurityCorrectionForUnfolding)
      {
          MakePhotonFigure30Like(
            hPhoRecoData,
            hPhoRecoSim,
            phoYieldDir,
            "leadingPhoton",
            "Leading-photon definition"
          );
      }

      TH2D* hPhoResp_measXtruth = TransposeTH2(
          hPhoRespSim,
          "h2_unfoldResponsePho_pTgamma_recoVsTruth",
          "Photon response matrix; p_{T}^{#gamma, reco} [GeV]; p_{T}^{#gamma, truth} [GeV]"
      );

      if (!hPhoResp_measXtruth)
      {
        cout << ANSI_BOLD_RED
             << "[ERROR] Failed to transpose photon response matrix. Aborting RooUnfold pipeline."
             << ANSI_RESET << "\n";
        delete hPhoRecoData;
        delete hPhoRecoSim;
        delete hPhoTruthSim;
        delete hPhoRespSim;
        return;
      }

        const int kDefaultBayesIterPho = 3;
        const int kMaxBayesIterPhoScan = 10;

        // Toy settings:
        //   - "final" is for the baseline unfolded spectrum that feeds your main outputs
        //   - "scan" is for iteration-stability / closure / half-closure diagnostics
        const int kNToysPhoFinal = 600;
        const int kNToysPhoScan  = 120;

        RooUnfoldResponse respPho(hPhoRecoSim, hPhoTruthSim, hPhoResp_measXtruth, "respPho", "respPho");

        const IterScanSummary phoIterScan = BuildPhotonIterScan(
          respPho,
          hPhoRecoData,
          hPhoTruthSim,
          kMaxBayesIterPhoScan,
          kNToysPhoScan
        );

        const int kBayesIterPho = (phoIterScan.bestIt > 0 ? phoIterScan.bestIt : kDefaultBayesIterPho);

        cout << ANSI_BOLD_CYN
             << "[PHO ITER AUTO] Using Bayes iteration " << kBayesIterPho;
        if (phoIterScan.bestIt > 0 && std::isfinite(phoIterScan.bestQuad))
        {
          cout << " from quadrature-sum minimum (" << phoIterScan.bestQuad << ")";
        }
        else
        {
          cout << " (fallback default; automatic scan unavailable)";
        }
        cout << ANSI_RESET << "\n";

        RooUnfoldBayes    unfoldPhoToy(&respPho, hPhoRecoData, kBayesIterPho);
        unfoldPhoToy.SetVerbose(0);
        unfoldPhoToy.SetNToys(kNToysPhoFinal);

        TH1* hPhoUnfoldTruth = nullptr;
        if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
        hPhoUnfoldTruth = unfoldPhoToy.Hreco(RooUnfold::kCovToy);
        if (gSystem) gSystem->RedirectOutput(0);
        if (hPhoUnfoldTruth) hPhoUnfoldTruth->SetDirectory(nullptr);

        RooUnfoldBayes    unfoldPhoCov(&respPho, hPhoRecoData, kBayesIterPho);
        unfoldPhoCov.SetVerbose(0);

        TH1* hPhoUnfoldTruth_cov = unfoldPhoCov.Hreco(RooUnfold::kCovariance);
        if (hPhoUnfoldTruth_cov) hPhoUnfoldTruth_cov->SetDirectory(nullptr);

        // Photon QA outputs
        {
          {
            TCanvas c("c_pho_resp","c_pho_resp", 900, 750);
            ApplyCanvasMargins2D(c);
            c.SetLogz();

            hPhoRespSim->SetTitle("");
            hPhoRespSim->GetXaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
            hPhoRespSim->GetYaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
            hPhoRespSim->Draw("colz");

            DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
            DrawLatexLines(0.14,0.78, { "SIM photon response", "truth #rightarrow reco" }, 0.030, 0.040);

            SaveCanvas(c, JoinPath(phoDir, "pho_response_truthVsReco.png"));
          }
          if (hPhoResp_measXtruth)
          {
              TCanvas c("c_pho_respT","c_pho_respT", 900, 750);
              ApplyCanvasMargins2D(c);
              c.SetLogz();

              hPhoResp_measXtruth->SetTitle("");
              hPhoResp_measXtruth->GetXaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
              hPhoResp_measXtruth->GetYaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");

              const auto& recoEdges  = kUnfoldRecoPtEdges;
              const auto& truthEdges = kUnfoldTruthPtEdges;

              const double recoMin     = (!recoEdges.empty() ? recoEdges.front() : 8.0);
              const double recoUFHi    = (recoEdges.size() >= 2 ? recoEdges[1] : 10.0);
              const double recoOFLo    = (recoEdges.size() >= 2 ? recoEdges[recoEdges.size() - 2] : 35.0);
              const double recoMax     = (!recoEdges.empty() ? recoEdges.back() : 40.0);
              const double recoDrawMax = recoOFLo;

              const double truthMin   = (!truthEdges.empty() ? truthEdges.front() : 5.0);
              const double truthUFMid = (truthEdges.size() >= 2 ? truthEdges[1] : 8.0);
              const double truthUFHi  = (truthEdges.size() >= 3 ? truthEdges[2] : 10.0);
              const double truthOFLo  = (truthEdges.size() >= 2 ? truthEdges[truthEdges.size() - 2] : 35.0);
              const double truthMax   = (!truthEdges.empty() ? truthEdges.back() : 40.0);

              // Show the full truth range, but cut the reco drawing range at 35 GeV
              // so the reco overflow support bin (35-40) is indicated by the dashed line
              // at 35 GeV rather than displayed as a full x-axis bin.
              hPhoResp_measXtruth->GetXaxis()->SetRangeUser(recoMin, recoDrawMax);
              hPhoResp_measXtruth->GetYaxis()->SetRangeUser(truthMin, truthMax);

              // Debug: verify the response histogram actually contains the expected bin edges
              cout << "  [pho_response_recoVsTruth] X(reco) nbins=" << hPhoResp_measXtruth->GetXaxis()->GetNbins()
                     << "  firstLowEdge=" << hPhoResp_measXtruth->GetXaxis()->GetBinLowEdge(1)
                     << "  lastUpEdge="   << hPhoResp_measXtruth->GetXaxis()->GetBinUpEdge(hPhoResp_measXtruth->GetXaxis()->GetNbins())
                     << "\n";
              cout << "  [pho_response_recoVsTruth] Y(truth) nbins=" << hPhoResp_measXtruth->GetYaxis()->GetNbins()
                     << "  firstLowEdge=" << hPhoResp_measXtruth->GetYaxis()->GetBinLowEdge(1)
                     << "  lastUpEdge="   << hPhoResp_measXtruth->GetYaxis()->GetBinUpEdge(hPhoResp_measXtruth->GetYaxis()->GetNbins())
                     << "\n";

              hPhoResp_measXtruth->GetZaxis()->SetTitle("Counts");
              hPhoResp_measXtruth->GetZaxis()->SetTitleSize(0.035);
              hPhoResp_measXtruth->GetZaxis()->SetTitleOffset(1.25);
              hPhoResp_measXtruth->GetZaxis()->SetLabelSize(0.030);

              hPhoResp_measXtruth->Draw("colz");
              if (gPad) gPad->Update();

              TLine lRecoUF(recoUFHi, truthMin, recoUFHi, truthMax);
              TLine lRecoOF(recoOFLo, truthMin, recoOFLo, truthMax);
              TLine lTruthUF1(recoMin, truthUFMid, recoDrawMax, truthUFMid);
              TLine lTruthUF2(recoMin, truthUFHi, recoDrawMax, truthUFHi);
              TLine lTruthOF(recoMin, truthOFLo, recoDrawMax, truthOFLo);

              lRecoUF.SetLineColor(kBlack);
              lRecoOF.SetLineColor(kBlack);
              lTruthUF1.SetLineColor(kBlack);
              lTruthUF2.SetLineColor(kBlack);
              lTruthOF.SetLineColor(kBlack);

              lRecoUF.SetLineStyle(2);
              lRecoOF.SetLineStyle(2);
              lTruthUF1.SetLineStyle(2);
              lTruthUF2.SetLineStyle(2);
              lTruthOF.SetLineStyle(2);

              lRecoUF.SetLineWidth(2);
              lRecoOF.SetLineWidth(2);
              lTruthUF1.SetLineWidth(2);
              lTruthUF2.SetLineWidth(2);
              lTruthOF.SetLineWidth(2);

              if (recoUFHi > recoMin && recoUFHi < recoMax) lRecoUF.Draw("same");
              if (recoOFLo > recoMin && recoOFLo < recoMax) lRecoOF.Draw("same");
              if (truthUFMid > truthMin && truthUFMid < truthMax) lTruthUF1.Draw("same");
              if (truthUFHi > truthMin && truthUFHi < truthMax) lTruthUF2.Draw("same");
              if (truthOFLo > truthMin && truthOFLo < truthMax) lTruthOF.Draw("same");

              DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
              DrawLatexLines(0.14,0.88, { "SIM photon response (transpose)", "reco #rightarrow truth axis order" }, 0.030, 0.040);
              DrawLatexLines(0.52,0.36,
                             {
                               "UF / OF support bins:",
                               "reco: 8-10 = UF, 35-40 = OF",
                               "truth: 5-8 and 8-10 = UF, 35-40 = OF"
                             },
                             0.026, 0.035);

              c.Modified();
              c.Update();

              SaveCanvas(c, JoinPath(phoDir, "pho_response_recoVsTruth.png"));
            }

          if (hPhoUnfoldTruth)
          {
            TH1* hRecoShape = CloneTH1(hPhoRecoData, "hPhoRecoData_forOverlay");
            TH1* hUnfShape  = CloneTH1(hPhoUnfoldTruth, "hPhoUnfoldTruth_forOverlay");
            if (hRecoShape && hUnfShape)
            {
              hRecoShape->SetDirectory(nullptr);
              hUnfShape->SetDirectory(nullptr);

              hRecoShape->SetLineColor(2);
              hRecoShape->SetMarkerColor(2);
              hRecoShape->SetMarkerStyle(20);
              hRecoShape->SetLineWidth(2);

              hUnfShape->SetLineColor(1);
              hUnfShape->SetMarkerColor(1);
              hUnfShape->SetMarkerStyle(24);
              hUnfShape->SetLineWidth(2);

              TCanvas c("c_pho_unf","c_pho_unf", 900, 700);
              ApplyCanvasMargins1D(c);

              const double maxv = std::max(hRecoShape->GetMaximum(), hUnfShape->GetMaximum());
              hRecoShape->SetMaximum(maxv * 1.35);
              hRecoShape->SetTitle("");
              hRecoShape->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
              hRecoShape->GetYaxis()->SetTitle("Counts");

              hRecoShape->Draw("E1");
              hUnfShape->Draw("E1 same");

              TLegend leg(0.55,0.76,0.92,0.90);
              leg.SetTextFont(42);
              leg.SetTextSize(0.032);
              leg.AddEntry(hRecoShape, "Run24pp Reco", "lep");
              leg.AddEntry(hUnfShape,  TString::Format("Unfolded (truth), Bayes it=%d", kBayesIterPho).Data(), "lep");
              leg.Draw();

              DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsData), 0.034, 0.045);
              DrawLatexLines(0.14,0.78, { "Photon unfolding: N_{#gamma}(p_{T}^{#gamma})" }, 0.030, 0.040);

              SaveCanvas(c, JoinPath(phoDir, "pho_unfolded_truth_pTgamma_overlay.png"));

              // Also save a log-y version (keep the linear-y output above as-is),
              // but only draw the common analysis bins (10-35 GeV) and add a ratio subpanel.
              {
                    const double xPlotMin = 10.0;
                    const double xPlotMax = 35.0;

                    auto scanMaxInRange = [&](TH1* h)->double
                    {
                      double maxY = 0.0;
                      if (!h) return maxY;

                      const int nb = h->GetNbinsX();
                      for (int ib = 1; ib <= nb; ++ib)
                      {
                        const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                        const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                        if (lo < xPlotMin || hi > xPlotMax) continue;

                        const double y  = h->GetBinContent(ib);
                        const double ey = h->GetBinError(ib);
                        if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                        if (y <= 0.0 && ey <= 0.0) continue;

                        const double v = y + ey;
                        if (v > maxY) maxY = v;
                      }
                      return maxY;
                    };

                    auto scanMinPosInRange = [&](TH1* h)->double
                    {
                      double minY = 1e99;
                      if (!h) return minY;

                      const int nb = h->GetNbinsX();
                      for (int ib = 1; ib <= nb; ++ib)
                      {
                        const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                        const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                        if (lo < xPlotMin || hi > xPlotMax) continue;

                        const double y = h->GetBinContent(ib);
                        if (!std::isfinite(y)) continue;
                        if (y > 0.0 && y < minY) minY = y;
                      }
                      return minY;
                    };

                    auto BuildPhotonAfterOverBeforeRatio = [&](TH1* hAfter, TH1* hBefore)->TH1D*
                    {
                      if (!hAfter || !hBefore) return nullptr;

                      vector<double> ratioEdges =
                      {
                        10.0, 12.0, 14.0, 16.0, 18.0,
                        20.0, 22.0, 24.0, 26.0, 35.0
                      };

                      TH1D* hRatio = new TH1D(
                        "hPho_afterOverBefore_ratio",
                        "",
                        (int)ratioEdges.size() - 1,
                        &ratioEdges[0]
                      );
                      hRatio->SetDirectory(nullptr);
                      hRatio->Sumw2();

                      for (int ib = 1; ib <= hRatio->GetNbinsX(); ++ib)
                      {
                        const double lo  = hRatio->GetXaxis()->GetBinLowEdge(ib);
                        const double hi  = hRatio->GetXaxis()->GetBinUpEdge(ib);
                        const double cen = hRatio->GetXaxis()->GetBinCenter(ib);

                        const int iAfter  = hAfter ->GetXaxis()->FindBin(cen);
                        const int iBefore = hBefore->GetXaxis()->FindBin(cen);

                        if (iAfter  < 1 || iAfter  > hAfter ->GetNbinsX()) continue;
                        if (iBefore < 1 || iBefore > hBefore->GetNbinsX()) continue;

                        const double afterLo  = hAfter ->GetXaxis()->GetBinLowEdge(iAfter);
                        const double afterHi  = hAfter ->GetXaxis()->GetBinUpEdge(iAfter);
                        const double beforeLo = hBefore->GetXaxis()->GetBinLowEdge(iBefore);
                        const double beforeHi = hBefore->GetXaxis()->GetBinUpEdge(iBefore);

                        if (std::fabs(afterLo  - lo) > 1e-6 || std::fabs(afterHi  - hi) > 1e-6) continue;
                        if (std::fabs(beforeLo - lo) > 1e-6 || std::fabs(beforeHi - hi) > 1e-6) continue;

                        const double num  = hAfter ->GetBinContent(iAfter);
                        const double eNum = hAfter ->GetBinError  (iAfter);
                        const double den  = hBefore->GetBinContent(iBefore);
                        const double eDen = hBefore->GetBinError  (iBefore);

                        if (!std::isfinite(num) || !std::isfinite(eNum) ||
                            !std::isfinite(den) || !std::isfinite(eDen) || den <= 0.0)
                        {
                          continue;
                        }

                        const double val = num / den;
                        const double var = (eNum * eNum) / (den * den)
                                         + (num * num * eDen * eDen) / (den * den * den * den);
                        const double err = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : 0.0;

                        hRatio->SetBinContent(ib, val);
                        hRatio->SetBinError  (ib, err);
                      }

                      return hRatio;
                    };

                    const double maxvTop = std::max(scanMaxInRange(hRecoShape), scanMaxInRange(hUnfShape));

                    double minPos = std::min(scanMinPosInRange(hRecoShape), scanMinPosInRange(hUnfShape));
                    if (!(minPos < 1e98)) minPos = 1e-3;

                    TH1D* hRatio = BuildPhotonAfterOverBeforeRatio(hUnfShape, hRecoShape);

                    double ratioMin = 0.8;
                    double ratioMax = 1.2;
                    if (hRatio)
                    {
                      double rMin =  1e99;
                      double rMax = -1e99;

                      for (int ib = 1; ib <= hRatio->GetNbinsX(); ++ib)
                      {
                        const double y  = hRatio->GetBinContent(ib);
                        const double ey = hRatio->GetBinError(ib);
                        if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                        if (y <= 0.0 && ey <= 0.0) continue;

                        rMin = std::min(rMin, y - ey);
                        rMax = std::max(rMax, y + ey);
                      }

                      if (rMin < 1e98 && rMax > -1e98 && rMin < rMax)
                      {
                        const double span = std::max(rMax - rMin, 0.04);
                        const double pad  = 0.18 * span;
                        ratioMin = rMin - pad;
                        ratioMax = rMax + pad;

                        if (ratioMin > 1.0) ratioMin = 1.0 - 0.5 * span;
                        if (ratioMax < 1.0) ratioMax = 1.0 + 0.5 * span;
                        if (ratioMin < 0.0) ratioMin = 0.0;
                      }
                    }

                    c.Clear();

                    TPad* pTop = new TPad("pTop_pho_unf_logy", "pTop_pho_unf_logy", 0.0, 0.30, 1.0, 1.0);
                    TPad* pBot = new TPad("pBot_pho_unf_logy", "pBot_pho_unf_logy", 0.0, 0.00, 1.0, 0.30);

                    pTop->SetLeftMargin(0.12);
                    pTop->SetRightMargin(0.04);
                    pTop->SetTopMargin(0.08);
                    pTop->SetBottomMargin(0.02);
                    pTop->SetTicks(1, 1);
                    pTop->SetLogy(1);

                    pBot->SetLeftMargin(0.12);
                    pBot->SetRightMargin(0.04);
                    pBot->SetTopMargin(0.02);
                    pBot->SetBottomMargin(0.30);
                    pBot->SetTicks(1, 1);
                    pBot->SetGridy(1);

                    pTop->Draw();
                    pBot->Draw();

                    pTop->cd();

                    hRecoShape->SetTitle("");
                    hRecoShape->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                    hUnfShape->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);

                    hRecoShape->GetXaxis()->SetTitle("");
                    hRecoShape->GetXaxis()->SetLabelSize(0.0);
                    hRecoShape->GetXaxis()->SetTitleSize(0.0);

                    hRecoShape->GetYaxis()->SetTitle("Counts");
                    hRecoShape->GetYaxis()->SetTitleSize(0.055);
                    hRecoShape->GetYaxis()->SetTitleOffset(0.95);
                    hRecoShape->GetYaxis()->SetLabelSize(0.045);

                    hRecoShape->SetMinimum(std::max(minPos * 0.5, 1e-6));
                    hRecoShape->SetMaximum((maxvTop > 0.0) ? (maxvTop * 3.0) : 1.0);

                    hRecoShape->Draw("E1");
                    hUnfShape->Draw("E1 same");

                    TLegend leg(0.55,0.76,0.92,0.90);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.032);
                    leg.AddEntry(hRecoShape, "PP DATA (reco)", "lep");
                    leg.AddEntry(hUnfShape,  TString::Format("Unfolded (truth), Bayes it=%d", kBayesIterPho).Data(), "lep");
                    leg.Draw();

                    DrawLatexLines(0.14,0.92,
                                   { "Photon 1D Unfolding, N_{#gamma}(p_{T}^{#gamma}), Photon 4 GeV + MBD NS #geq 1" },
                                   0.034, 0.045);

                    pBot->cd();

                    if (hRatio)
                    {
                      hRatio->SetTitle("");
                      hRatio->SetMarkerStyle(20);
                      hRatio->SetMarkerSize(0.95);
                      hRatio->SetLineWidth(2);
                      hRatio->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                      hRatio->GetYaxis()->SetTitle("After / before unfolding");
                      hRatio->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                      hRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
                      hRatio->GetXaxis()->SetTitleSize(0.12);
                      hRatio->GetXaxis()->SetLabelSize(0.11);
                      hRatio->GetXaxis()->SetTitleOffset(1.00);
                      hRatio->GetYaxis()->SetTitleSize(0.10);
                      hRatio->GetYaxis()->SetLabelSize(0.09);
                      hRatio->GetYaxis()->SetTitleOffset(0.55);
                      hRatio->GetYaxis()->SetNdivisions(505);
                      hRatio->Draw("E1");

                      TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");
                    }
                    else
                    {
                      TH1F frame("frame_pho_ratio","", 1, xPlotMin, xPlotMax);
                      frame.SetMinimum(ratioMin);
                      frame.SetMaximum(ratioMax);
                      frame.SetTitle("");
                      frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                      frame.GetYaxis()->SetTitle("After / before unfolding");
                      frame.GetXaxis()->SetTitleSize(0.12);
                      frame.GetXaxis()->SetLabelSize(0.11);
                      frame.GetXaxis()->SetTitleOffset(1.00);
                      frame.GetYaxis()->SetTitleSize(0.10);
                      frame.GetYaxis()->SetLabelSize(0.09);
                      frame.GetYaxis()->SetTitleOffset(0.55);
                      frame.GetYaxis()->SetNdivisions(505);
                      frame.Draw("axis");

                      TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");
                    }

                    c.cd();
                    c.Modified();
                    c.Update();
                    SaveCanvas(c, JoinPath(phoDir, "pho_unfolded_truth_pTgamma_overlay_logy.png"));

                  if (hRatio) delete hRatio;
              }

              // -----------------------------------------------------------------
              // 3-curve overlay: raw  vs  purity-corrected  vs  unfolded (logy)
              // Only meaningful when purity correction was applied so that we have
              // distinct raw and corrected inputs.
              // Output: <phoDir>/pho_raw_purityCorr_unfolded_pTgamma_overlay_logy.png
              // -----------------------------------------------------------------
              if (gApplyPurityCorrectionForUnfolding)
              {
                  const double xPlotMin = 10.0;
                  const double xPlotMax = 35.0;

                  TH1* hRaw = CloneTH1(hPhoRecoData_in,   "hPhoRecoRaw_3ov");
                  TH1* hPur = CloneTH1(hPhoRecoData,       "hPhoRecoPur_3ov");
                  TH1* hUnf = CloneTH1(hPhoUnfoldTruth,    "hPhoUnfold_3ov");

                  if (hRaw && hPur && hUnf)
                  {
                    hRaw->SetDirectory(nullptr);
                    hPur->SetDirectory(nullptr);
                    hUnf->SetDirectory(nullptr);

                    hRaw->SetLineColor(kGray + 2);
                    hRaw->SetMarkerColor(kGray + 2);
                    hRaw->SetMarkerStyle(21);
                    hRaw->SetMarkerSize(0.9);
                    hRaw->SetLineWidth(2);

                    hPur->SetLineColor(kRed);
                    hPur->SetMarkerColor(kRed);
                    hPur->SetMarkerStyle(20);
                    hPur->SetMarkerSize(0.9);
                    hPur->SetLineWidth(2);

                    hUnf->SetLineColor(kBlue + 1);
                    hUnf->SetMarkerColor(kBlue + 1);
                    hUnf->SetMarkerStyle(24);
                    hUnf->SetMarkerSize(0.9);
                    hUnf->SetLineWidth(2);

                    auto scanMax3 = [&](TH1* h)->double
                    {
                      double mx = 0.0;
                      if (!h) return mx;
                      for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                      {
                        const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                        const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                        if (lo < xPlotMin || hi > xPlotMax) continue;
                        const double y  = h->GetBinContent(ib);
                        const double ey = h->GetBinError(ib);
                        if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                        if (y <= 0.0 && ey <= 0.0) continue;
                        mx = std::max(mx, y + ey);
                      }
                      return mx;
                    };

                    auto scanMinPos3 = [&](TH1* h)->double
                    {
                      double mn = 1e99;
                      if (!h) return mn;
                      for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                      {
                        const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                        const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                        if (lo < xPlotMin || hi > xPlotMax) continue;
                        const double y = h->GetBinContent(ib);
                        if (std::isfinite(y) && y > 0.0 && y < mn) mn = y;
                      }
                      return mn;
                    };

                    const double maxvTop3 = std::max({scanMax3(hRaw), scanMax3(hPur), scanMax3(hUnf)});
                    double minPos3 = std::min({scanMinPos3(hRaw), scanMinPos3(hPur), scanMinPos3(hUnf)});
                    if (!(minPos3 < 1e98)) minPos3 = 1e-3;

                    // -- build ratio hists (denominator = raw) --
                    auto buildRatio3 = [&](TH1* hNum, TH1* hDen, const char* name)->TH1D*
                    {
                      if (!hNum || !hDen) return nullptr;
                      vector<double> edges = { 10.0, 12.0, 14.0, 16.0, 18.0,
                                               20.0, 22.0, 24.0, 26.0, 35.0 };
                      TH1D* hR = new TH1D(name, "", (int)edges.size() - 1, &edges[0]);
                      hR->SetDirectory(nullptr);
                      hR->Sumw2();
                      for (int ib = 1; ib <= hR->GetNbinsX(); ++ib)
                      {
                        const double cen = hR->GetXaxis()->GetBinCenter(ib);
                        const int iN = hNum->GetXaxis()->FindBin(cen);
                        const int iD = hDen->GetXaxis()->FindBin(cen);
                        if (iN < 1 || iN > hNum->GetNbinsX()) continue;
                        if (iD < 1 || iD > hDen->GetNbinsX()) continue;
                        const double num  = hNum->GetBinContent(iN);
                        const double eNum = hNum->GetBinError(iN);
                        const double den  = hDen->GetBinContent(iD);
                        const double eDen = hDen->GetBinError(iD);
                        if (!std::isfinite(num) || !std::isfinite(den) || den <= 0.0) continue;
                        const double val = num / den;
                        const double var = (eNum*eNum)/(den*den)
                                         + (num*num*eDen*eDen)/(den*den*den*den);
                        hR->SetBinContent(ib, val);
                        hR->SetBinError(ib, (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : 0.0);
                      }
                      return hR;
                    };

                    TH1D* hRatioPurOverRaw = buildRatio3(hPur, hRaw, "hRatio3_purOverRaw");
                    TH1D* hRatioUnfOverRaw = buildRatio3(hUnf, hRaw, "hRatio3_unfOverRaw");

                    // -- auto-range ratio panel --
                    double ratioMin3 = 0.5, ratioMax3 = 1.5;
                    {
                      double rLo =  1e99, rHi = -1e99;
                      auto scanR = [&](TH1D* h)
                      {
                        if (!h) return;
                        for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                        {
                          const double y  = h->GetBinContent(ib);
                          const double ey = h->GetBinError(ib);
                          if (!std::isfinite(y) || (y <= 0.0 && ey <= 0.0)) continue;
                          rLo = std::min(rLo, y - ey);
                          rHi = std::max(rHi, y + ey);
                        }
                      };
                      scanR(hRatioPurOverRaw);
                      scanR(hRatioUnfOverRaw);
                      if (rLo < 1e98 && rHi > -1e98 && rLo < rHi)
                      {
                        const double span = std::max(rHi - rLo, 0.04);
                        const double pad  = 0.18 * span;
                        ratioMin3 = rLo - pad;
                        ratioMax3 = rHi + pad;
                        if (ratioMin3 > 1.0) ratioMin3 = 1.0 - 0.5 * span;
                        if (ratioMax3 < 1.0) ratioMax3 = 1.0 + 0.5 * span;
                        if (ratioMin3 < 0.0) ratioMin3 = 0.0;
                      }
                    }

                    TCanvas c3("c_pho_3ov_logy", "c_pho_3ov_logy", 900, 800);

                    TPad* pTop3 = new TPad("pTop_3ov", "pTop_3ov", 0.0, 0.30, 1.0, 1.0);
                    TPad* pBot3 = new TPad("pBot_3ov", "pBot_3ov", 0.0, 0.00, 1.0, 0.30);

                    pTop3->SetLeftMargin(0.12);
                    pTop3->SetRightMargin(0.04);
                    pTop3->SetTopMargin(0.08);
                    pTop3->SetBottomMargin(0.02);
                    pTop3->SetTicks(1, 1);
                    pTop3->SetLogy(1);

                    pBot3->SetLeftMargin(0.12);
                    pBot3->SetRightMargin(0.04);
                    pBot3->SetTopMargin(0.02);
                    pBot3->SetBottomMargin(0.30);
                    pBot3->SetTicks(1, 1);
                    pBot3->SetGridy(1);

                    pTop3->Draw();
                    pBot3->Draw();

                    // --- top: logy overlay of 3 curves ---
                    pTop3->cd();

                    hRaw->SetTitle("");
                    hRaw->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                    hPur->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                    hUnf->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);

                    hRaw->GetXaxis()->SetTitle("");
                    hRaw->GetXaxis()->SetLabelSize(0.0);
                    hRaw->GetXaxis()->SetTitleSize(0.0);

                    hRaw->GetYaxis()->SetTitle("Counts");
                    hRaw->GetYaxis()->SetTitleSize(0.055);
                    hRaw->GetYaxis()->SetTitleOffset(0.95);
                    hRaw->GetYaxis()->SetLabelSize(0.045);

                    hRaw->SetMinimum(std::max(minPos3 * 0.5, 1e-6));
                    hRaw->SetMaximum((maxvTop3 > 0.0) ? (maxvTop3 * 3.0) : 1.0);

                    hRaw->Draw("E1");
                    hPur->Draw("E1 same");
                    hUnf->Draw("E1 same");

                    TLegend leg3(0.48, 0.70, 0.92, 0.90);
                    leg3.SetTextFont(42);
                    leg3.SetTextSize(0.032);
                    leg3.AddEntry(hRaw, "PP DATA (reco, non-purity-corrected)", "lep");
                    leg3.AddEntry(hPur, "PP DATA (reco, purity-corrected)", "lep");
                    leg3.AddEntry(hUnf, TString::Format("Unfolded (truth), Bayes it=%d", kBayesIterPho).Data(), "lep");
                    leg3.Draw();

                    DrawLatexLines(0.14, 0.92,
                                   { "Photon 1D: raw vs purity-corrected vs unfolded, N_{#gamma}(p_{T}^{#gamma})" },
                                   0.034, 0.045);

                    // --- bottom: ratios over raw ---
                    pBot3->cd();

                    bool drewFirst = false;
                    if (hRatioPurOverRaw)
                    {
                      hRatioPurOverRaw->SetTitle("");
                      hRatioPurOverRaw->SetMarkerStyle(20);
                      hRatioPurOverRaw->SetMarkerSize(0.85);
                      hRatioPurOverRaw->SetMarkerColor(kRed);
                      hRatioPurOverRaw->SetLineColor(kRed);
                      hRatioPurOverRaw->SetLineWidth(2);
                      hRatioPurOverRaw->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                      hRatioPurOverRaw->GetYaxis()->SetTitle("Ratio / raw");
                      hRatioPurOverRaw->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                      hRatioPurOverRaw->GetYaxis()->SetRangeUser(ratioMin3, ratioMax3);
                      hRatioPurOverRaw->GetXaxis()->SetTitleSize(0.12);
                      hRatioPurOverRaw->GetXaxis()->SetLabelSize(0.11);
                      hRatioPurOverRaw->GetXaxis()->SetTitleOffset(1.00);
                      hRatioPurOverRaw->GetYaxis()->SetTitleSize(0.10);
                      hRatioPurOverRaw->GetYaxis()->SetLabelSize(0.09);
                      hRatioPurOverRaw->GetYaxis()->SetTitleOffset(0.55);
                      hRatioPurOverRaw->GetYaxis()->SetNdivisions(505);
                      hRatioPurOverRaw->Draw("E1");
                      drewFirst = true;
                    }

                    if (hRatioUnfOverRaw)
                    {
                      hRatioUnfOverRaw->SetMarkerStyle(24);
                      hRatioUnfOverRaw->SetMarkerSize(0.85);
                      hRatioUnfOverRaw->SetMarkerColor(kBlue + 1);
                      hRatioUnfOverRaw->SetLineColor(kBlue + 1);
                      hRatioUnfOverRaw->SetLineWidth(2);

                      if (!drewFirst)
                      {
                        hRatioUnfOverRaw->SetTitle("");
                        hRatioUnfOverRaw->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hRatioUnfOverRaw->GetYaxis()->SetTitle("Ratio / raw");
                        hRatioUnfOverRaw->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                        hRatioUnfOverRaw->GetYaxis()->SetRangeUser(ratioMin3, ratioMax3);
                        hRatioUnfOverRaw->GetXaxis()->SetTitleSize(0.12);
                        hRatioUnfOverRaw->GetXaxis()->SetLabelSize(0.11);
                        hRatioUnfOverRaw->GetXaxis()->SetTitleOffset(1.00);
                        hRatioUnfOverRaw->GetYaxis()->SetTitleSize(0.10);
                        hRatioUnfOverRaw->GetYaxis()->SetLabelSize(0.09);
                        hRatioUnfOverRaw->GetYaxis()->SetTitleOffset(0.55);
                        hRatioUnfOverRaw->GetYaxis()->SetNdivisions(505);
                        hRatioUnfOverRaw->Draw("E1");
                      }
                      else
                      {
                        hRatioUnfOverRaw->Draw("E1 same");
                      }
                    }

                    if (!drewFirst && !hRatioUnfOverRaw)
                    {
                      TH1F frame3("frame_3ov_ratio","", 1, xPlotMin, xPlotMax);
                      frame3.SetMinimum(ratioMin3);
                      frame3.SetMaximum(ratioMax3);
                      frame3.SetTitle("");
                      frame3.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                      frame3.GetYaxis()->SetTitle("Ratio / raw");
                      frame3.GetXaxis()->SetTitleSize(0.12);
                      frame3.GetXaxis()->SetLabelSize(0.11);
                      frame3.GetXaxis()->SetTitleOffset(1.00);
                      frame3.GetYaxis()->SetTitleSize(0.10);
                      frame3.GetYaxis()->SetLabelSize(0.09);
                      frame3.GetYaxis()->SetTitleOffset(0.55);
                      frame3.GetYaxis()->SetNdivisions(505);
                      frame3.Draw("axis");
                    }

                    {
                      TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");
                    }

                    // ratio legend
                    {
                      TLegend legR(0.48, 0.72, 0.92, 0.92);
                      legR.SetBorderSize(0);
                      legR.SetFillStyle(0);
                      legR.SetTextFont(42);
                      legR.SetTextSize(0.085);
                      if (hRatioPurOverRaw) legR.AddEntry(hRatioPurOverRaw, "purity-corr / raw", "lep");
                      if (hRatioUnfOverRaw) legR.AddEntry(hRatioUnfOverRaw, "unfolded / raw",    "lep");
                      legR.DrawClone();
                    }

                    c3.cd();
                    c3.Modified();
                    c3.Update();
                    SaveCanvas(c3, JoinPath(phoDir, "pho_raw_purityCorr_unfolded_pTgamma_overlay_logy.png"));

                    if (hRatioPurOverRaw) delete hRatioPurOverRaw;
                    if (hRatioUnfOverRaw) delete hRatioUnfOverRaw;
                  }

                  if (hRaw) delete hRaw;
                  if (hPur) delete hPur;
                  if (hUnf) delete hUnf;
                }

                // photon efficiency/purity diagnostics
                //
                // These are SIM truth-matching bookkeeping diagnostics that explain the
                // normalization gap in before/after unfolding overlays.
                //
                // (1) Purity (SIM, reco space):
                //     purity(pT) = 1 - N_fakeReco(pT)/N_reco(pT)
                //
                // (2) Efficiency (SIM, truth space):
                //     eff(pT)    = 1 - N_missTruth(pT)/N_truth(pT)
                //
                // (3) Truth/Reconstruction scale factor (SIM):
                //     N_truth(pT)/N_reco(pT)  (with bin-mapping when axes differ)
                //
                // (4) Implied "efficiency-like" curve in DATA (not purely data-driven):
                //     eps_eff,data(pT) = N_reco,data(pT) / N_truth,data(unfolded)(pT)
                //     (computed with bin-mapping when axes differ)
                //
                // Outputs (to <phoDir>, i.e. unfolding/radii/<rXX>/photons):
                //   - pho_efficiencyEff_data_vs_pTgamma.png
                //   - pho_purity_sim_vs_pTgamma.png                 (if fakes hist exists)
                //   - pho_efficiency_sim_vs_pTgamma.png             (if misses hist exists)
                //   - pho_truthOverReco_sim_vs_pTgamma.png          (bin-mapped)
                //   - pho_recoOverTruth_sim_vs_pTgamma.png          (bin-mapped)
                //   - pho_efficiencyEff_data_vs_efficiency_sim.png  (overlay, if both exist)
                // -------------------------------------------------------------------
                {
                  auto sameBinning = [&](TH1* a, TH1* b)->bool
                  {
                    if (!a || !b) return false;
                    if (a->GetNbinsX() != b->GetNbinsX()) return false;
                    const int nb = a->GetNbinsX();
                    for (int ib = 1; ib <= nb + 1; ++ib)
                    {
                      const double ea = a->GetXaxis()->GetBinUpEdge(ib);
                      const double eb = b->GetXaxis()->GetBinUpEdge(ib);
                      if (std::fabs(ea - eb) > 1e-9) return false;
                    }
                    return true;
                  };

                  auto mapToRefBinning = [&](TH1* src, TH1* ref, const char* newName)->TH1*
                  {
                    if (!src || !ref) return nullptr;

                    TH1* h = CloneTH1(ref, newName);
                    if (!h) return nullptr;

                    h->SetDirectory(nullptr);
                    EnsureSumw2(h);
                    h->Reset("ICES");

                    const int nb = ref->GetNbinsX();
                    int nBad = 0;

                    for (int ib = 1; ib <= nb; ++ib)
                    {
                      const double x  = ref->GetXaxis()->GetBinCenter(ib);
                      const int isrc  = src->GetXaxis()->FindBin(x);

                      const double xsLo = src->GetXaxis()->GetBinLowEdge(isrc);
                      const double xsHi = src->GetXaxis()->GetBinUpEdge(isrc);
                      const double xrLo = ref->GetXaxis()->GetBinLowEdge(ib);
                      const double xrHi = ref->GetXaxis()->GetBinUpEdge(ib);

                      if (std::fabs(xsLo - xrLo) > 1e-3 || std::fabs(xsHi - xrHi) > 1e-3) ++nBad;

                      h->SetBinContent(ib, src->GetBinContent(isrc));
                      h->SetBinError  (ib, src->GetBinError  (isrc));
                    }

                    if (nBad > 0)
                    {
                      cout << ANSI_BOLD_YEL
                           << "[WARN] Photon diagnostics: mapped histogram '" << newName
                           << "' had " << nBad << " bins with mismatched edges (copied by bin-center)."
                           << ANSI_RESET << "\n";
                    }

                    return h;
                  };

                  auto drawLineAtOne = [&](TH1* h)
                  {
                    if (!h) return;
                    const double xmin = h->GetXaxis()->GetXmin();
                    const double xmax = h->GetXaxis()->GetXmax();
                    TLine l1(xmin, 1.0, xmax, 1.0);
                    l1.SetLineStyle(2);
                    l1.SetLineWidth(2);
                    l1.Draw("same");
                  };

                  // -----------------------------
                  // DATA: eps_eff = N_reco,data / N_truth,data(unfolded)
                  // -----------------------------
                  TH1* hEffData = nullptr;
                  {
                    TH1* hRecoD  = hPhoRecoData;
                    TH1* hTruthD = hPhoUnfoldTruth;

                    if (hRecoD && hTruthD)
                    {
                      TH1* hRecoD_m = nullptr;

                      if (sameBinning(hRecoD, hTruthD))
                      {
                        hRecoD_m = CloneTH1(hRecoD, "h_pho_recoData_counts_forEff");
                        if (hRecoD_m) { hRecoD_m->SetDirectory(nullptr); EnsureSumw2(hRecoD_m); }
                      }
                      else
                      {
                        hRecoD_m = mapToRefBinning(hRecoD, hTruthD, "h_pho_recoData_counts_mappedToTruthBins_forEff");
                      }

                      if (hRecoD_m)
                      {
                        hEffData = CloneTH1(hRecoD_m, "h_pho_effData_recoOverUnfoldTruth");
                        if (hEffData)
                        {
                          hEffData->SetDirectory(nullptr);
                          EnsureSumw2(hEffData);
                          hEffData->Divide(hTruthD);

                          hEffData->SetTitle("");
                          hEffData->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                          hEffData->GetYaxis()->SetTitle("#epsilon_{#gamma}^{eff,data} = N_{#gamma}^{reco,data} / N_{#gamma}^{truth,data (unfolded)}");
                          hEffData->SetMarkerStyle(20);
                          hEffData->SetMarkerSize(1.1);
                          hEffData->SetLineWidth(2);

                          TCanvas c("c_pho_effData", "c_pho_effData", 900, 700);
                          ApplyCanvasMargins1D(c);

                          hEffData->GetYaxis()->SetRangeUser(0.0, 1.2);
                          hEffData->Draw("E1");
                          drawLineAtOne(hEffData);

                          DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsData), 0.034, 0.045);
                          DrawLatexLines(0.14,0.78, { "Photon unfolding diagnostic: implied #epsilon_{#gamma}^{eff,data}(p_{T}^{#gamma})" }, 0.030, 0.040);

                          SaveCanvas(c, JoinPath(phoDir, "pho_efficiencyEff_data_vs_pTgamma.png"));
                        }
                      }

                      if (hRecoD_m && hRecoD_m != hRecoD) delete hRecoD_m;
                    }
                  }

                  // -----------------------------
                  // SIM: purity, efficiency, truth/reco ratios
                  // -----------------------------
                  TH1* hPurSim = nullptr;
                  TH1* hEffSim = nullptr;
                  TH1* hTruthOverRecoSim = nullptr;
                  TH1* hRecoOverTruthSim = nullptr;

                    // Purity: 1 - fakes/reco  (reco space)
                    if (hPhoRecoSim && hPhoRecoFakesSim_in)
                    {
                      TH1* hFakeOverReco = CloneTH1(hPhoRecoFakesSim_in, "h_pho_fakeOverReco_sim");
                      if (hFakeOverReco)
                      {
                        hFakeOverReco->SetDirectory(nullptr);
                        EnsureSumw2(hFakeOverReco);
                        hFakeOverReco->Divide(hPhoRecoSim);

                      hPurSim = CloneTH1(hFakeOverReco, "h_pho_purity_sim");
                      if (hPurSim)
                      {
                        hPurSim->SetDirectory(nullptr);
                        EnsureSumw2(hPurSim);

                        for (int ib = 1; ib <= hPurSim->GetNbinsX(); ++ib)
                        {
                          const double v  = 1.0 - hFakeOverReco->GetBinContent(ib);
                          const double ev = hFakeOverReco->GetBinError(ib);
                          hPurSim->SetBinContent(ib, v);
                          hPurSim->SetBinError  (ib, ev);
                        }

                        hPurSim->SetTitle("");
                        hPurSim->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hPurSim->GetYaxis()->SetTitle("Purity(p_{T}^{#gamma}) = 1 - N_{fake}^{reco}/N_{#gamma}^{reco}");
                        hPurSim->SetMarkerStyle(21);
                        hPurSim->SetMarkerSize(1.1);
                        hPurSim->SetLineWidth(2);

                        TCanvas c("c_pho_purity_sim", "c_pho_purity_sim", 900, 700);
                        ApplyCanvasMargins1D(c);

                        hPurSim->GetYaxis()->SetRangeUser(0.0, 1.2);
                        hPurSim->Draw("E1");
                        drawLineAtOne(hPurSim);

                        DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
                        DrawLatexLines(0.14,0.78, { "SIM photon purity: 1 - (reco fakes)/(reco selected)" }, 0.030, 0.040);

                        SaveCanvas(c, JoinPath(phoDir, "pho_purity_sim_vs_pTgamma.png"));
                      }

                      delete hFakeOverReco;
                    }
                  }
                  else
                  {
                    cout << ANSI_BOLD_YEL
                         << "[WARN] Photon purity plot: missing SIM fakes histogram (h_unfoldRecoPhoFakes_pTgamma). Skipping purity plot."
                         << ANSI_RESET << "\n";
                  }

                    // Efficiency: 1 - misses/truth (truth space)
                  if (hPhoTruthSim && hPhoTruthMissesSim_in)
                  {
                      TH1* hMissOverTruth = CloneTH1(hPhoTruthMissesSim_in, "h_pho_missOverTruth_sim");
                      if (hMissOverTruth)
                      {
                        hMissOverTruth->SetDirectory(nullptr);
                        EnsureSumw2(hMissOverTruth);
                        hMissOverTruth->Divide(hPhoTruthSim);

                      hEffSim = CloneTH1(hMissOverTruth, "h_pho_efficiency_sim");
                      if (hEffSim)
                      {
                        hEffSim->SetDirectory(nullptr);
                        EnsureSumw2(hEffSim);

                        for (int ib = 1; ib <= hEffSim->GetNbinsX(); ++ib)
                        {
                          const double v  = 1.0 - hMissOverTruth->GetBinContent(ib);
                          const double ev = hMissOverTruth->GetBinError(ib);
                          hEffSim->SetBinContent(ib, v);
                          hEffSim->SetBinError  (ib, ev);
                        }

                        hEffSim->SetTitle("");
                        hEffSim->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hEffSim->GetYaxis()->SetTitle("#epsilon_{#gamma}^{MC}(p_{T}^{#gamma}) = 1 - N_{miss}^{truth}/N_{#gamma}^{truth}");
                        hEffSim->SetMarkerStyle(22);
                        hEffSim->SetMarkerSize(1.1);
                        hEffSim->SetLineWidth(2);

                        TCanvas c("c_pho_efficiency_sim", "c_pho_efficiency_sim", 900, 700);
                        ApplyCanvasMargins1D(c);

                        hEffSim->GetYaxis()->SetRangeUser(0.0, 1.2);
                        hEffSim->Draw("E1");
                        drawLineAtOne(hEffSim);

                        DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
                        DrawLatexLines(0.14,0.78, { "SIM photon efficiency: 1 - (truth misses)/(truth signal)" }, 0.030, 0.040);

                        SaveCanvas(c, JoinPath(phoDir, "pho_efficiency_sim_vs_pTgamma.png"));
                      }

                      delete hMissOverTruth;
                    }
                  }
                  else
                  {
                    cout << ANSI_BOLD_YEL
                         << "[WARN] Photon efficiency plot: missing SIM misses histogram (h_unfoldTruthPhoMisses_pTgamma). Skipping efficiency plot."
                         << ANSI_RESET << "\n";
                  }

                  // Truth/Reco ratios in SIM (bin-mapped so it is well-defined)
                  if (hPhoTruthSim && hPhoRecoSim)
                  {
                    TH1* hRecoSim_mTruth = nullptr;
                    if (sameBinning(hPhoRecoSim, hPhoTruthSim))
                    {
                      hRecoSim_mTruth = CloneTH1(hPhoRecoSim, "h_pho_recoSim_counts_forTruthOverReco");
                      if (hRecoSim_mTruth) { hRecoSim_mTruth->SetDirectory(nullptr); EnsureSumw2(hRecoSim_mTruth); }
                    }
                    else
                    {
                      hRecoSim_mTruth = mapToRefBinning(hPhoRecoSim, hPhoTruthSim, "h_pho_recoSim_counts_mappedToTruthBins_forTruthOverReco");
                    }

                    if (hRecoSim_mTruth)
                    {
                      hTruthOverRecoSim = CloneTH1(hPhoTruthSim, "h_pho_truthOverReco_sim");
                      if (hTruthOverRecoSim)
                      {
                        hTruthOverRecoSim->SetDirectory(nullptr);
                        EnsureSumw2(hTruthOverRecoSim);
                        hTruthOverRecoSim->Divide(hRecoSim_mTruth);

                        hTruthOverRecoSim->SetTitle("");
                        hTruthOverRecoSim->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hTruthOverRecoSim->GetYaxis()->SetTitle("N_{#gamma}^{truth} / N_{#gamma}^{reco}  (SIM, bin-mapped)");
                        hTruthOverRecoSim->SetMarkerStyle(20);
                        hTruthOverRecoSim->SetMarkerSize(1.1);
                        hTruthOverRecoSim->SetLineWidth(2);

                        TCanvas c("c_pho_truthOverReco_sim", "c_pho_truthOverReco_sim", 900, 700);
                        ApplyCanvasMargins1D(c);

                        hTruthOverRecoSim->Draw("E1");

                        DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
                        DrawLatexLines(0.14,0.78, { "SIM normalization lever-arm: N_{#gamma}^{truth} / N_{#gamma}^{reco}" }, 0.030, 0.040);

                        SaveCanvas(c, JoinPath(phoDir, "pho_truthOverReco_sim_vs_pTgamma.png"));
                      }

                      if (hRecoSim_mTruth && hRecoSim_mTruth != hPhoRecoSim) delete hRecoSim_mTruth;
                    }

                    TH1* hTruthSim_mReco = nullptr;
                    if (sameBinning(hPhoTruthSim, hPhoRecoSim))
                    {
                      hTruthSim_mReco = CloneTH1(hPhoTruthSim, "h_pho_truthSim_counts_forRecoOverTruth");
                      if (hTruthSim_mReco) { hTruthSim_mReco->SetDirectory(nullptr); EnsureSumw2(hTruthSim_mReco); }
                    }
                    else
                    {
                      hTruthSim_mReco = mapToRefBinning(hPhoTruthSim, hPhoRecoSim, "h_pho_truthSim_counts_mappedToRecoBins_forRecoOverTruth");
                    }

                    if (hTruthSim_mReco)
                    {
                      hRecoOverTruthSim = CloneTH1(hPhoRecoSim, "h_pho_recoOverTruth_sim");
                      if (hRecoOverTruthSim)
                      {
                        hRecoOverTruthSim->SetDirectory(nullptr);
                        EnsureSumw2(hRecoOverTruthSim);
                        hRecoOverTruthSim->Divide(hTruthSim_mReco);

                        hRecoOverTruthSim->SetTitle("");
                        hRecoOverTruthSim->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hRecoOverTruthSim->GetYaxis()->SetTitle("N_{#gamma}^{reco} / N_{#gamma}^{truth}  (SIM, bin-mapped)");
                        hRecoOverTruthSim->SetMarkerStyle(20);
                        hRecoOverTruthSim->SetMarkerSize(1.1);
                        hRecoOverTruthSim->SetLineWidth(2);

                        TCanvas c("c_pho_recoOverTruth_sim", "c_pho_recoOverTruth_sim", 900, 700);
                        ApplyCanvasMargins1D(c);

                        hRecoOverTruthSim->GetYaxis()->SetRangeUser(0.0, 1.2);
                        hRecoOverTruthSim->Draw("E1");
                        drawLineAtOne(hRecoOverTruthSim);

                        DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
                        DrawLatexLines(0.14,0.78, { "SIM efficiency-like: N_{#gamma}^{reco} / N_{#gamma}^{truth}" }, 0.030, 0.040);

                        SaveCanvas(c, JoinPath(phoDir, "pho_recoOverTruth_sim_vs_pTgamma.png"));
                      }

                      if (hTruthSim_mReco && hTruthSim_mReco != hPhoTruthSim) delete hTruthSim_mReco;
                    }
                  }

                  // -----------------------------
                  // Overlay: eps_eff,data vs eps_MC (if both exist)
                  // -----------------------------
                  if (hEffData && hEffSim)
                  {
                    TCanvas c("c_pho_effData_vs_effSim", "c_pho_effData_vs_effSim", 900, 700);
                    ApplyCanvasMargins1D(c);

                    TH1* hFrame = CloneTH1(hEffSim, "hFrame_effData_vs_effSim");
                    if (hFrame)
                    {
                      hFrame->SetDirectory(nullptr);
                      hFrame->Reset("ICES");
                      hFrame->SetTitle("");
                      hFrame->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                      hFrame->GetYaxis()->SetTitle("#epsilon_{#gamma}(p_{T}^{#gamma})");
                      hFrame->GetYaxis()->SetRangeUser(0.0, 1.2);
                      hFrame->Draw("axis");
                    }

                    hEffSim->SetMarkerColor(kBlue + 1);
                    hEffSim->SetLineColor(kBlue + 1);
                    hEffData->SetMarkerColor(kBlack);
                    hEffData->SetLineColor(kBlack);

                    hEffSim->Draw("E1 same");
                    hEffData->Draw("E1 same");

                    TGraphErrors gLegData(1), gLegSim(1);
                    {
                      int ib = 1;
                      const double x  = hEffData->GetXaxis()->GetBinCenter(ib);
                      const double y  = hEffData->GetBinContent(ib);
                      const double ey = hEffData->GetBinError(ib);
                      gLegData.SetPoint(0, x, y);
                      gLegData.SetPointError(0, 0.0, ey);
                      gLegData.SetMarkerStyle(hEffData->GetMarkerStyle());
                      gLegData.SetMarkerSize(hEffData->GetMarkerSize());
                      gLegData.SetMarkerColor(hEffData->GetMarkerColor());
                      gLegData.SetLineColor(hEffData->GetLineColor());
                      gLegData.SetLineWidth(hEffData->GetLineWidth());
                    }
                    {
                      int ib = 1;
                      const double x  = hEffSim->GetXaxis()->GetBinCenter(ib);
                      const double y  = hEffSim->GetBinContent(ib);
                      const double ey = hEffSim->GetBinError(ib);
                      gLegSim.SetPoint(0, x, y);
                      gLegSim.SetPointError(0, 0.0, ey);
                      gLegSim.SetMarkerStyle(hEffSim->GetMarkerStyle());
                      gLegSim.SetMarkerSize(hEffSim->GetMarkerSize());
                      gLegSim.SetMarkerColor(hEffSim->GetMarkerColor());
                      gLegSim.SetLineColor(hEffSim->GetLineColor());
                      gLegSim.SetLineWidth(hEffSim->GetLineWidth());
                    }

                      TLegend leg(0.22, 0.74, 0.45, 0.88);
                      leg.SetBorderSize(0);
                      leg.SetFillStyle(0);
                      leg.SetTextFont(42);
                      leg.SetTextSize(0.035);
                      leg.AddEntry(&gLegData, "#epsilon_{#gamma}^{eff,data} (reco / unfolded truth)", "pe");
                      leg.AddEntry(&gLegSim,  "#epsilon_{#gamma}^{MC} (1 - misses/truth)", "pe");
                      leg.Draw();

                      {
                        const double xmin = hFrame ? hFrame->GetXaxis()->GetXmin() : hEffSim->GetXaxis()->GetXmin();
                        const double xmax = hFrame ? hFrame->GetXaxis()->GetXmax() : hEffSim->GetXaxis()->GetXmax();
                        TLine l1(xmin, 1.0, xmax, 1.0);
                        l1.SetLineStyle(2);
                        l1.SetLineWidth(2);
                        l1.SetLineColor(kRed + 1);
                        l1.Draw("same");
                      }

                      {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(22);
                        tx.SetTextSize(0.040);
                        tx.DrawLatex(0.50, 0.965, "Photon efficiency diagnostic: DATA reco / unfolded truth vs SIM (1 - misses / truth)");
                    }

                    SaveCanvas(c, JoinPath(phoDir, "pho_efficiencyEff_data_vs_efficiency_sim.png"));

                    if (hFrame) delete hFrame;
                  }

                  if (hRecoOverTruthSim) delete hRecoOverTruthSim;
                  if (hTruthOverRecoSim) delete hTruthOverRecoSim;
                  if (hEffSim) delete hEffSim;
                  if (hPurSim) delete hPurSim;
                  if (hEffData) delete hEffData;
                }

                delete hRecoShape;
                delete hUnfShape;
            }
          }
        }
        // Photon QA outputs
        {
          auto sameBinning = [&](TH1* a, TH1* b)->bool
          {
            if (!a || !b) return false;
            if (a->GetNbinsX() != b->GetNbinsX()) return false;
            const int nb = a->GetNbinsX();
            for (int ib = 1; ib <= nb + 1; ++ib)
            {
              const double ea = a->GetXaxis()->GetBinUpEdge(ib);
              const double eb = b->GetXaxis()->GetBinUpEdge(ib);
              if (std::fabs(ea - eb) > 1e-9) return false;
            }
            return true;
          };

          auto mapToRefBinning = [&](TH1* src, TH1* ref, const char* newName)->TH1*
          {
            if (!src || !ref) return nullptr;

            TH1* h = CloneTH1(ref, newName);
            if (!h) return nullptr;

            h->SetDirectory(nullptr);
            EnsureSumw2(h);
            h->Reset("ICES");

            const int nb = ref->GetNbinsX();
            int nBad = 0;

            for (int ib = 1; ib <= nb; ++ib)
            {
              const double x  = ref->GetXaxis()->GetBinCenter(ib);
              const int isrc  = src->GetXaxis()->FindBin(x);

              const double xsLo = src->GetXaxis()->GetBinLowEdge(isrc);
              const double xsHi = src->GetXaxis()->GetBinUpEdge(isrc);
              const double xrLo = ref->GetXaxis()->GetBinLowEdge(ib);
              const double xrHi = ref->GetXaxis()->GetBinUpEdge(ib);

              if (std::fabs(xsLo - xrLo) > 1e-3 || std::fabs(xsHi - xrHi) > 1e-3) ++nBad;

              h->SetBinContent(ib, src->GetBinContent(isrc));
              h->SetBinError  (ib, src->GetBinError  (isrc));
            }

            if (nBad > 0)
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Photon QA: mapped histogram '" << newName
                   << "' had " << nBad << " bins with mismatched edges (copied by bin-center)."
                   << ANSI_RESET << "\n";
            }

            return h;
          };

          auto cloneOrMapToRef = [&](TH1* src, TH1* ref, const char* newName)->TH1*
          {
            if (!src || !ref) return nullptr;

            TH1* h = nullptr;
            if (sameBinning(src, ref))
            {
              h = CloneTH1(src, newName);
              if (h)
              {
                h->SetDirectory(nullptr);
                EnsureSumw2(h);
              }
            }
            else
            {
              h = mapToRefBinning(src, ref, newName);
            }
            return h;
          };

          auto makeRatioOnRef = [&](TH1* hNum, TH1* hDen, TH1* hRef, const char* newName)->TH1*
          {
            if (!hNum || !hDen || !hRef) return nullptr;

            TH1* hNumR = cloneOrMapToRef(hNum, hRef, TString::Format("%s_numRef", newName).Data());
            TH1* hDenR = cloneOrMapToRef(hDen, hRef, TString::Format("%s_denRef", newName).Data());
            if (!hNumR || !hDenR)
            {
              if (hNumR) delete hNumR;
              if (hDenR) delete hDenR;
              return nullptr;
            }

            TH1* hR = CloneTH1(hNumR, newName);
            if (hR)
            {
              hR->SetDirectory(nullptr);
              EnsureSumw2(hR);
              hR->Divide(hDenR);
            }

            delete hNumR;
            delete hDenR;
            return hR;
          };

          auto makeOneMinus = [&](TH1* hIn, const char* newName)->TH1*
          {
            if (!hIn) return nullptr;

            TH1* h = CloneTH1(hIn, newName);
            if (!h) return nullptr;

            h->SetDirectory(nullptr);
            EnsureSumw2(h);

            for (int ib = 0; ib <= h->GetNbinsX() + 1; ++ib)
            {
              if (ib == 0 || ib == h->GetNbinsX() + 1)
              {
                h->SetBinContent(ib, 0.0);
                h->SetBinError  (ib, 0.0);
                continue;
              }

              const double v  = 1.0 - hIn->GetBinContent(ib);
              const double ev = hIn->GetBinError(ib);
              h->SetBinContent(ib, v);
              h->SetBinError  (ib, ev);
            }

            return h;
          };

          auto makeHistFromAxis = [&](const TAxis* ax, const char* newName, const char* yTitle)->TH1D*
          {
            if (!ax) return nullptr;

            TH1D* h = nullptr;
            if (ax->GetXbins() && ax->GetXbins()->GetSize() > 0)
            {
              h = new TH1D(newName, "", ax->GetNbins(), ax->GetXbins()->GetArray());
            }
            else
            {
              h = new TH1D(newName, "", ax->GetNbins(), ax->GetXmin(), ax->GetXmax());
            }

            if (!h) return nullptr;
            h->SetDirectory(nullptr);
            h->Sumw2();
            h->SetTitle("");
            h->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
            h->GetYaxis()->SetTitle(yTitle);
            return h;
          };

          auto normalizeResponseByTruth = [&](TH2* hIn, const char* newName)->TH2*
          {
            if (!hIn) return nullptr;

            TH2* h = CloneTH2(hIn, newName);
            if (!h) return nullptr;

            h->SetDirectory(nullptr);
            EnsureSumw2(h);
            h->Reset("ICES");

            for (int ix = 1; ix <= hIn->GetNbinsX(); ++ix)
            {
              double sum = 0.0;
              for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy) sum += hIn->GetBinContent(ix, iy);
              if (!(sum > 0.0)) continue;

              for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy)
              {
                const double v = hIn->GetBinContent(ix, iy) / sum;
                h->SetBinContent(ix, iy, v);
                h->SetBinError  (ix, iy, 0.0);
              }
            }

            return h;
          };

          auto normalizeResponseByReco = [&](TH2* hIn, const char* newName)->TH2*
          {
            if (!hIn) return nullptr;

            TH2* h = CloneTH2(hIn, newName);
            if (!h) return nullptr;

            h->SetDirectory(nullptr);
            EnsureSumw2(h);
            h->Reset("ICES");

            for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy)
            {
              double sum = 0.0;
              for (int ix = 1; ix <= hIn->GetNbinsX(); ++ix) sum += hIn->GetBinContent(ix, iy);
              if (!(sum > 0.0)) continue;

              for (int ix = 1; ix <= hIn->GetNbinsX(); ++ix)
              {
                const double v = hIn->GetBinContent(ix, iy) / sum;
                h->SetBinContent(ix, iy, v);
                h->SetBinError  (ix, iy, 0.0);
              }
            }

            return h;
          };

          auto normalizeRecoFeedIn = [&](TH2* hIn, const char* newName)->TH2*
          {
            if (!hIn) return nullptr;

            TH2* h = CloneTH2(hIn, newName);
            if (!h) return nullptr;

            h->SetDirectory(nullptr);
            EnsureSumw2(h);
            h->Reset("ICES");

            for (int ix = 1; ix <= hIn->GetNbinsX(); ++ix)
            {
              double sum = 0.0;
              for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy) sum += hIn->GetBinContent(ix, iy);
              if (!(sum > 0.0)) continue;

              for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy)
              {
                const double v = hIn->GetBinContent(ix, iy) / sum;
                h->SetBinContent(ix, iy, v);
                h->SetBinError  (ix, iy, 0.0);
              }
            }

            return h;
          };

          auto buildResponseMeanRecoOverTruth = [&](TH2* hResp, const char* newName)->TH1D*
          {
            if (!hResp) return nullptr;

            TH1D* h = makeHistFromAxis(hResp->GetXaxis(), newName, "Mean(reco / truth)");
            if (!h) return nullptr;

            for (int ix = 1; ix <= hResp->GetNbinsX(); ++ix)
            {
              const double xTruth = hResp->GetXaxis()->GetBinCenter(ix);
              if (!(xTruth > 0.0)) continue;

              double sumW = 0.0;
              double sumR = 0.0;
              double sumR2 = 0.0;

              for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy)
              {
                const double w = hResp->GetBinContent(ix, iy);
                if (!(w > 0.0)) continue;

                const double xReco = hResp->GetYaxis()->GetBinCenter(iy);
                const double r = xReco / xTruth;

                sumW  += w;
                sumR  += w * r;
                sumR2 += w * r * r;
              }

              if (!(sumW > 0.0)) continue;

              const double mean = sumR / sumW;
              const double var  = std::max(0.0, sumR2 / sumW - mean * mean);
              const double err  = std::sqrt(var / sumW);

              h->SetBinContent(ix, mean);
              h->SetBinError  (ix, err);
            }

            return h;
          };

          auto buildResponseWidthRecoOverTruth = [&](TH2* hResp, const char* newName)->TH1D*
          {
            if (!hResp) return nullptr;

            TH1D* h = makeHistFromAxis(hResp->GetXaxis(), newName, "RMS(reco / truth)");
            if (!h) return nullptr;

            for (int ix = 1; ix <= hResp->GetNbinsX(); ++ix)
            {
              const double xTruth = hResp->GetXaxis()->GetBinCenter(ix);
              if (!(xTruth > 0.0)) continue;

              double sumW = 0.0;
              double sumR = 0.0;
              double sumR2 = 0.0;

              for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy)
              {
                const double w = hResp->GetBinContent(ix, iy);
                if (!(w > 0.0)) continue;

                const double xReco = hResp->GetYaxis()->GetBinCenter(iy);
                const double r = xReco / xTruth;

                sumW  += w;
                sumR  += w * r;
                sumR2 += w * r * r;
              }

              if (!(sumW > 0.0)) continue;

              const double mean = sumR / sumW;
              const double var  = std::max(0.0, sumR2 / sumW - mean * mean);
              const double rms  = std::sqrt(var);
              const double err  = (sumW > 1.0) ? (rms / std::sqrt(2.0 * (sumW - 1.0))) : 0.0;

              h->SetBinContent(ix, rms);
              h->SetBinError  (ix, err);
            }

            return h;
          };

          auto buildResponseDiagonalFraction = [&](TH2* hResp, const char* newName)->TH1D*
          {
            if (!hResp) return nullptr;

            TH1D* h = makeHistFromAxis(hResp->GetXaxis(), newName, "Same-bin fraction");
            if (!h) return nullptr;

            for (int ix = 1; ix <= hResp->GetNbinsX(); ++ix)
            {
              const double xLo = hResp->GetXaxis()->GetBinLowEdge(ix);
              const double xHi = hResp->GetXaxis()->GetBinUpEdge(ix);

              int iyMatch = -1;
              for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy)
              {
                const double yLo = hResp->GetYaxis()->GetBinLowEdge(iy);
                const double yHi = hResp->GetYaxis()->GetBinUpEdge(iy);
                if (std::fabs(yLo - xLo) < 1e-9 && std::fabs(yHi - xHi) < 1e-9)
                {
                  iyMatch = iy;
                  break;
                }
              }

              double sum = 0.0;
              for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy) sum += hResp->GetBinContent(ix, iy);
              if (!(sum > 0.0) || iyMatch < 0) continue;

              const double frac = hResp->GetBinContent(ix, iyMatch) / sum;
              const double err  = std::sqrt(std::max(0.0, frac * (1.0 - frac) / sum));

              h->SetBinContent(ix, frac);
              h->SetBinError  (ix, err);
            }

            return h;
          };

          auto findExactBinByEdges = [&](const TAxis* ax, double lo, double hi)->int
          {
            if (!ax) return -1;
            for (int ib = 1; ib <= ax->GetNbins(); ++ib)
            {
              const double bLo = ax->GetBinLowEdge(ib);
              const double bHi = ax->GetBinUpEdge(ib);
              if (std::fabs(bLo - lo) < 1e-9 && std::fabs(bHi - hi) < 1e-9) return ib;
            }
            return -1;
          };

          auto styleHist = [&](TH1* h, int color, int marker)->void
          {
            if (!h) return;
            h->SetLineColor(color);
            h->SetMarkerColor(color);
            h->SetMarkerStyle(marker);
            h->SetMarkerSize(0.95);
            h->SetLineWidth(2);
            h->SetTitle("");
          };

          auto prepPad = [&]()->void
          {
            if (!gPad) return;
            gPad->SetLeftMargin(0.12);
            gPad->SetRightMargin(0.05);
            gPad->SetBottomMargin(0.12);
            gPad->SetTopMargin(0.10);
            gPad->SetTicks(1, 1);
          };

          auto drawPadTitle = [&](const string& txt)->void
          {
            TLatex tx;
            tx.SetNDC();
            tx.SetTextFont(42);
            tx.SetTextAlign(22);
            tx.SetTextSize(0.050);
            tx.DrawLatex(0.50, 0.95, txt.c_str());
          };

          auto drawCanvasTitle = [&](const string& txt, double size)->void
          {
            TLatex tx;
            tx.SetNDC();
            tx.SetTextFont(42);
            tx.SetTextAlign(22);
            tx.SetTextSize(size);
            tx.DrawLatex(0.50, 0.985, txt.c_str());
          };

          auto drawCanvasNote = [&](const vector<string>& lines, double x, double y, int align, double size)->void
          {
            TLatex tx;
            tx.SetNDC();
            tx.SetTextFont(42);
            tx.SetTextAlign(align);
            tx.SetTextSize(size);
            double yy = y;
            for (const auto& line : lines)
            {
              tx.DrawLatex(x, yy, line.c_str());
              yy -= 0.045;
            }
          };

          auto drawMissingPad = [&](const string& yTitle, const string& msg)->void
          {
            TH1F frame("frame_missing_phoQA", "", 1, 10.0, 35.0);
            frame.SetMinimum(0.0);
            frame.SetMaximum(1.0);
            frame.SetTitle("");
            frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
            frame.GetYaxis()->SetTitle(yTitle.c_str());
            frame.Draw("axis");

            TLatex tx;
            tx.SetNDC();
            tx.SetTextFont(42);
            tx.SetTextAlign(22);
            tx.SetTextSize(0.050);
            tx.DrawLatex(0.50, 0.50, msg.c_str());
          };

          auto setRatioRange = [&](TH1* h, double yMinDefault, double yMaxDefault)->void
          {
            if (!h) return;

            double yMin =  1e99;
            double yMax = -1e99;
            bool have = false;

            for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
            {
              const double lo = h->GetXaxis()->GetBinLowEdge(ib);
              const double hi = h->GetXaxis()->GetBinUpEdge(ib);
              if (lo < 10.0 || hi > 35.0) continue;

              const double y  = h->GetBinContent(ib);
              const double ey = h->GetBinError(ib);
              if (!std::isfinite(y) || !std::isfinite(ey)) continue;
              if (y <= 0.0 && ey <= 0.0) continue;

              have = true;
              yMin = std::min(yMin, y - ey);
              yMax = std::max(yMax, y + ey);
            }

            if (!have || !(yMin < yMax))
            {
              h->GetYaxis()->SetRangeUser(yMinDefault, yMaxDefault);
              return;
            }

            const double pad = std::max(0.08, 0.18 * (yMax - yMin));
            yMin = std::min(yMin, 1.0) - pad;
            yMax = std::max(yMax, 1.0) + pad;
            if (yMin < 0.0) yMin = 0.0;
            if (yMax <= yMin) yMax = yMin + 0.5;
            h->GetYaxis()->SetRangeUser(yMin, yMax);
          };

          const bool buildPhotonQA = gApplyPurityCorrectionForUnfolding;
          const string phoQADir       = JoinPath(phoDir, "QA");
          const string qaMasterDir    = JoinPath(phoQADir, "01_master");
          const string qaBudgetDir    = JoinPath(phoQADir, "02_budget");
          const string qaResponseDir  = JoinPath(phoQADir, "03_response");
          const string qaSupportDir   = JoinPath(phoQADir, "04_support");

          if (buildPhotonQA)
          {
            EnsureDir(phoQADir);
            EnsureDir(qaMasterDir);
            EnsureDir(qaBudgetDir);
            EnsureDir(qaResponseDir);
            EnsureDir(qaSupportDir);
          }

          if (hPhoResp_measXtruth)
          {
            TCanvas c("c_pho_respT","c_pho_respT", 900, 750);
            ApplyCanvasMargins2D(c);
            c.SetLogz();

            hPhoResp_measXtruth->SetTitle("");
            hPhoResp_measXtruth->GetXaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
            hPhoResp_measXtruth->GetXaxis()->SetTitleOffset(1.15);
            hPhoResp_measXtruth->GetYaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
            hPhoResp_measXtruth->GetZaxis()->SetTitle("Counts");

            const auto& recoEdges  = kUnfoldRecoPtEdges;
            const auto& truthEdges = kUnfoldTruthPtEdges;

            const double recoMin     = (!recoEdges.empty() ? recoEdges.front() : 8.0);
            const double recoUFHi    = (recoEdges.size() >= 2 ? recoEdges[1] : 10.0);
            const double recoOFLo    = (recoEdges.size() >= 2 ? recoEdges[recoEdges.size() - 2] : 35.0);
            const double recoDrawMax = recoOFLo;

            const double truthMin   = (!truthEdges.empty() ? truthEdges.front() : 5.0);
            const double truthUFMid = (truthEdges.size() >= 2 ? truthEdges[1] : 8.0);
            const double truthUFHi  = (truthEdges.size() >= 3 ? truthEdges[2] : 10.0);
            const double truthOFLo  = (truthEdges.size() >= 2 ? truthEdges[truthEdges.size() - 2] : 35.0);
            const double truthMax   = (!truthEdges.empty() ? truthEdges.back() : 40.0);

            hPhoResp_measXtruth->GetXaxis()->SetRangeUser(recoMin, recoDrawMax);
            hPhoResp_measXtruth->GetYaxis()->SetRangeUser(truthMin, truthMax);
            hPhoResp_measXtruth->Draw("colz");
            if (gPad) gPad->Update();

            TLine lRecoUF(recoUFHi, truthMin, recoUFHi, truthMax);
            TLine lRecoOF(recoOFLo, truthMin, recoOFLo, truthMax);
            TLine lTruthUF1(recoMin, truthUFMid, recoDrawMax, truthUFMid);
            TLine lTruthUF2(recoMin, truthUFHi, recoDrawMax, truthUFHi);
            TLine lTruthOF(recoMin, truthOFLo, recoDrawMax, truthOFLo);

            lRecoUF.SetLineColor(kBlack);
            lRecoOF.SetLineColor(kBlack);
            lTruthUF1.SetLineColor(kBlack);
            lTruthUF2.SetLineColor(kBlack);
            lTruthOF.SetLineColor(kBlack);

            lRecoUF.SetLineStyle(2);
            lRecoOF.SetLineStyle(2);
            lTruthUF1.SetLineStyle(2);
            lTruthUF2.SetLineStyle(2);
            lTruthOF.SetLineStyle(2);

            lRecoUF.SetLineWidth(2);
            lRecoOF.SetLineWidth(2);
            lTruthUF1.SetLineWidth(2);
            lTruthUF2.SetLineWidth(2);
            lTruthOF.SetLineWidth(2);

            if (recoUFHi > recoMin) lRecoUF.Draw("same");
            if (recoOFLo > recoMin) lRecoOF.Draw("same");
            if (truthUFMid > truthMin) lTruthUF1.Draw("same");
            if (truthUFHi > truthMin) lTruthUF2.Draw("same");
            if (truthOFLo > truthMin) lTruthOF.Draw("same");

            DrawLatexLines(0.14,0.95, DefaultHeaderLines(dsSim), 0.034, 0.045);
            DrawLatexLines(0.52,0.36,
                           {
                             "UF / OF support bins:",
                             "reco: 8-10 = UF, 35-40 = OF",
                             "truth: 5-8 and 8-10 = UF, 35-40 = OF"
                           },
                           0.026, 0.035);

            c.Modified();
            c.Update();
            SaveCanvas(c, JoinPath(phoDir, "pho_response_recoVsTruth.png"));
          }

          if (hPhoUnfoldTruth)
          {
            TH1* hRecoShape = CloneTH1(hPhoRecoData, "hPhoRecoData_forOverlayLog");
            TH1* hUnfShape  = CloneTH1(hPhoUnfoldTruth, "hPhoUnfoldTruth_forOverlayLog");

            if (hRecoShape && hUnfShape)
            {
              hRecoShape->SetDirectory(nullptr);
              hUnfShape->SetDirectory(nullptr);
              EnsureSumw2(hRecoShape);
              EnsureSumw2(hUnfShape);

              hRecoShape->SetLineColor(2);
              hRecoShape->SetMarkerColor(2);
              hRecoShape->SetMarkerStyle(20);
              hRecoShape->SetLineWidth(2);

              hUnfShape->SetLineColor(1);
              hUnfShape->SetMarkerColor(1);
              hUnfShape->SetMarkerStyle(24);
              hUnfShape->SetLineWidth(2);

              const double xPlotMin = 10.0;
              const double xPlotMax = 35.0;

              auto scanMaxInRange = [&](TH1* h)->double
              {
                double maxY = 0.0;
                if (!h) return maxY;

                const int nb = h->GetNbinsX();
                for (int ib = 1; ib <= nb; ++ib)
                {
                  const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                  const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                  if (lo < xPlotMin || hi > xPlotMax) continue;

                  const double y  = h->GetBinContent(ib);
                  const double ey = h->GetBinError(ib);
                  if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                  if (y <= 0.0 && ey <= 0.0) continue;

                  const double v = y + ey;
                  if (v > maxY) maxY = v;
                }
                return maxY;
              };

              auto scanMinPosInRange = [&](TH1* h)->double
              {
                double minY = 1e99;
                if (!h) return minY;

                const int nb = h->GetNbinsX();
                for (int ib = 1; ib <= nb; ++ib)
                {
                  const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                  const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                  if (lo < xPlotMin || hi > xPlotMax) continue;

                  const double y = h->GetBinContent(ib);
                  if (!std::isfinite(y)) continue;
                  if (y > 0.0 && y < minY) minY = y;
                }
                return minY;
              };

              const double maxvTop = std::max(scanMaxInRange(hRecoShape), scanMaxInRange(hUnfShape));
              double minPos = std::min(scanMinPosInRange(hRecoShape), scanMinPosInRange(hUnfShape));
              if (!(minPos < 1e98)) minPos = 1e-3;

              TH1* hRatio = makeRatioOnRef(hPhoUnfoldTruth, hPhoRecoData, hPhoRecoData, "hPho_afterOverBefore_ratio");
              if (hRatio)
              {
                hRatio->SetDirectory(nullptr);
                EnsureSumw2(hRatio);
              }

              TCanvas c("c_pho_unf_logy","c_pho_unf_logy", 900, 700);
              TPad* pTop = new TPad("pTop_pho_unf_logy", "pTop_pho_unf_logy", 0.0, 0.30, 1.0, 1.0);
              TPad* pBot = new TPad("pBot_pho_unf_logy", "pBot_pho_unf_logy", 0.0, 0.00, 1.0, 0.30);

              pTop->SetLeftMargin(0.12);
              pTop->SetRightMargin(0.04);
              pTop->SetTopMargin(0.08);
              pTop->SetBottomMargin(0.02);
              pTop->SetTicks(1, 1);
              pTop->SetLogy(1);

              pBot->SetLeftMargin(0.12);
              pBot->SetRightMargin(0.04);
              pBot->SetTopMargin(0.02);
              pBot->SetBottomMargin(0.30);
              pBot->SetTicks(1, 1);
              pBot->SetGridy(1);

              pTop->Draw();
              pBot->Draw();

              pTop->cd();
              hRecoShape->SetTitle("");
              hRecoShape->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
              hUnfShape->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
              hRecoShape->GetXaxis()->SetTitle("");
              hRecoShape->GetXaxis()->SetLabelSize(0.0);
              hRecoShape->GetXaxis()->SetTitleSize(0.0);
              hRecoShape->GetYaxis()->SetTitle("Counts");
              hRecoShape->GetYaxis()->SetTitleSize(0.055);
              hRecoShape->GetYaxis()->SetTitleOffset(0.95);
              hRecoShape->GetYaxis()->SetLabelSize(0.045);
              hRecoShape->SetMinimum(std::max(minPos * 0.5, 1e-6));
              hRecoShape->SetMaximum((maxvTop > 0.0) ? (maxvTop * 3.0) : 1.0);
              hRecoShape->Draw("E1");
              hUnfShape->Draw("E1 same");

              TLegend leg(0.55,0.76,0.92,0.90);
              leg.SetTextFont(42);
              leg.SetTextSize(0.032);
              leg.AddEntry(hRecoShape, "PP DATA (reco)", "lep");
              leg.AddEntry(hUnfShape,  TString::Format("Unfolded (truth), Bayes it=%d", kBayesIterPho).Data(), "lep");
              leg.Draw();

              DrawLatexLines(0.14,0.92,
                             { "Photon 1D Unfolding, N_{#gamma}(p_{T}^{#gamma}), Photon 4 GeV + MBD NS #geq 1" },
                             0.034, 0.045);

              pBot->cd();
              if (hRatio)
              {
                hRatio->SetTitle("");
                hRatio->SetMarkerStyle(20);
                hRatio->SetMarkerSize(0.95);
                hRatio->SetLineWidth(2);
                hRatio->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                hRatio->GetYaxis()->SetTitle("After / before unfolding");
                hRatio->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                setRatioRange(hRatio, 0.8, 1.2);
                hRatio->GetXaxis()->SetTitleSize(0.12);
                hRatio->GetXaxis()->SetLabelSize(0.11);
                hRatio->GetXaxis()->SetTitleOffset(1.00);
                hRatio->GetYaxis()->SetTitleSize(0.10);
                hRatio->GetYaxis()->SetLabelSize(0.09);
                hRatio->GetYaxis()->SetTitleOffset(0.55);
                hRatio->GetYaxis()->SetNdivisions(505);
                hRatio->Draw("E1");

                TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                l1.SetLineStyle(2);
                l1.SetLineWidth(2);
                l1.Draw("same");
              }
              else
              {
                TH1F frame("frame_pho_ratio","", 1, xPlotMin, xPlotMax);
                frame.SetMinimum(0.8);
                frame.SetMaximum(1.2);
                frame.SetTitle("");
                frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                frame.GetYaxis()->SetTitle("After / before unfolding");
                frame.GetXaxis()->SetTitleSize(0.12);
                frame.GetXaxis()->SetLabelSize(0.11);
                frame.GetXaxis()->SetTitleOffset(1.00);
                frame.GetYaxis()->SetTitleSize(0.10);
                frame.GetYaxis()->SetLabelSize(0.09);
                frame.GetYaxis()->SetTitleOffset(0.55);
                frame.GetYaxis()->SetNdivisions(505);
                frame.Draw("axis");

                TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                l1.SetLineStyle(2);
                l1.SetLineWidth(2);
                l1.Draw("same");
              }

              c.cd();
              c.Modified();
              c.Update();
              SaveCanvas(c, JoinPath(phoDir, "pho_unfolded_truth_pTgamma_overlay_logy.png"));

              if (hRatio) delete hRatio;
              delete pTop;
              delete pBot;
            }

            if (hRecoShape) delete hRecoShape;
            if (hUnfShape) delete hUnfShape;
          }

          if (buildPhotonQA)
          {
            TH1* hRecoDataQA   = cloneOrMapToRef(hPhoRecoData, hPhoRecoData, "hPhoRecoData_QA");
            TH1* hRecoSimQA    = cloneOrMapToRef(hPhoRecoSim,  hPhoRecoData, "hPhoRecoSim_QA");
            TH1* hTruthSimQA   = cloneOrMapToRef(hPhoTruthSim, hPhoRecoData, "hPhoTruthSim_onReco_QA");
            TH1* hMissesSimQA  = cloneOrMapToRef(hPhoTruthMissesSim_in, hPhoRecoData, "hPhoTruthMisses_onReco_QA");
            TH1* hUnfoldDataQA = cloneOrMapToRef(hPhoUnfoldTruth, hPhoRecoData, "hPhoUnfoldTruth_onReco_QA");

            TH1* hTruthOverRecoSim = makeRatioOnRef(hPhoTruthSim, hPhoRecoSim, hPhoRecoData, "hPho_truthOverRecoSim_QA");
            TH1* hRecoOverData     = makeRatioOnRef(hPhoRecoSim, hPhoRecoData, hPhoRecoData, "hPho_recoSimOverData_QA");
            TH1* hTruthOverData    = makeRatioOnRef(hPhoTruthSim, hPhoRecoData, hPhoRecoData, "hPho_truthSimOverData_QA");

            TH1* hFakeFracSim = makeRatioOnRef(hPhoRecoFakesSim_in, hPhoRecoSim, hPhoRecoData, "hPho_fakeFractionSim_QA");
            TH1* hPuritySim   = makeOneMinus(hFakeFracSim, "hPho_puritySim_QA");
            TH1* hMissFracSim = makeRatioOnRef(hPhoTruthMissesSim_in, hPhoTruthSim, hPhoRecoData, "hPho_missFractionSim_QA");
            TH1* hEffSim      = makeOneMinus(hMissFracSim, "hPho_efficiencySim_QA");
            TH1* hAfterOverBefore = makeRatioOnRef(hPhoUnfoldTruth, hPhoRecoData, hPhoRecoData, "hPho_afterOverBefore_QA");

            styleHist(hRecoDataQA,   kBlack,    20);
            styleHist(hRecoSimQA,    kBlue + 1, 24);
            styleHist(hTruthSimQA,   kRed + 1,  25);
            styleHist(hTruthOverRecoSim, kBlack,    20);
            styleHist(hRecoOverData,     kBlue + 1, 20);
            styleHist(hTruthOverData,    kRed + 1,  20);
            styleHist(hFakeFracSim,      kRed + 1,  20);
            styleHist(hPuritySim,        kBlue + 1, 20);
            styleHist(hMissFracSim,      kRed + 1,  20);
            styleHist(hEffSim,           kBlue + 1, 20);
            styleHist(hAfterOverBefore,  kBlack,    20);

            if (hRecoDataQA && hRecoSimQA && hTruthSimQA && hTruthOverRecoSim && hRecoOverData && hTruthOverData)
            {
              TCanvas c("c_pho_masterQA", "c_pho_masterQA", 1500, 1050);
              c.Divide(2, 2, 0.001, 0.001);

              c.cd(1);
              prepPad();
              if (gPad) gPad->SetLogy(1);

              double maxTop = 0.0;
              auto scanMasterMax = [&](TH1* h)->void
              {
                if (!h) return;
                for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                {
                  const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                  const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                  if (lo < 10.0 || hi > 35.0) continue;
                  const double v = h->GetBinContent(ib) + h->GetBinError(ib);
                  if (std::isfinite(v) && v > maxTop) maxTop = v;
                }
              };
              scanMasterMax(hRecoDataQA);
              scanMasterMax(hRecoSimQA);
              scanMasterMax(hTruthSimQA);

              double minPos = 1e99;
              auto scanMasterMinPos = [&](TH1* h)->void
              {
                if (!h) return;
                for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                {
                  const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                  const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                  if (lo < 10.0 || hi > 35.0) continue;
                  const double v = h->GetBinContent(ib);
                  if (std::isfinite(v) && v > 0.0 && v < minPos) minPos = v;
                }
              };
              scanMasterMinPos(hRecoDataQA);
              scanMasterMinPos(hRecoSimQA);
              scanMasterMinPos(hTruthSimQA);
              if (!(minPos < 1e98)) minPos = 1e-3;

              hRecoDataQA->SetMinimum(std::max(0.5 * minPos, 1e-6));
              hRecoDataQA->SetMaximum((maxTop > 0.0) ? (3.0 * maxTop) : 1.0);
              hRecoDataQA->GetXaxis()->SetRangeUser(10.0, 35.0);
              hRecoDataQA->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
              hRecoDataQA->GetYaxis()->SetTitle("Counts");
              hRecoDataQA->Draw("E1");
              hRecoSimQA->GetXaxis()->SetRangeUser(10.0, 35.0);
              hRecoSimQA->Draw("E1 same");
              hTruthSimQA->GetXaxis()->SetRangeUser(10.0, 35.0);
              hTruthSimQA->Draw("E1 same");
              drawPadTitle("DATA reco vs SIM reco vs SIM truth");

              TLegend leg(0.52, 0.70, 0.90, 0.88);
              leg.SetBorderSize(0);
              leg.SetFillStyle(0);
              leg.SetTextFont(42);
              leg.SetTextSize(0.040);
              leg.AddEntry(hRecoDataQA, "DATA reco", "pe");
              leg.AddEntry(hRecoSimQA,  "SIM reco",  "pe");
              leg.AddEntry(hTruthSimQA, "SIM truth (bin-mapped)", "pe");
              leg.Draw();

              c.cd(2);
              prepPad();
              hTruthOverRecoSim->GetXaxis()->SetRangeUser(10.0, 35.0);
              hTruthOverRecoSim->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
              hTruthOverRecoSim->GetYaxis()->SetTitle("SIM truth / SIM reco");
              setRatioRange(hTruthOverRecoSim, 0.5, 2.5);
              hTruthOverRecoSim->Draw("E1");
              {
                TLine l1(10.0, 1.0, 35.0, 1.0);
                l1.SetLineStyle(2);
                l1.SetLineWidth(2);
                l1.Draw("same");
              }
              drawPadTitle("SIM lever-arm");

              c.cd(3);
              prepPad();
              hRecoOverData->GetXaxis()->SetRangeUser(10.0, 35.0);
              hRecoOverData->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
              hRecoOverData->GetYaxis()->SetTitle("SIM reco / DATA reco");
              setRatioRange(hRecoOverData, 0.5, 1.5);
              hRecoOverData->Draw("E1");
              {
                TLine l1(10.0, 1.0, 35.0, 1.0);
                l1.SetLineStyle(2);
                l1.SetLineWidth(2);
                l1.Draw("same");
              }
              drawPadTitle("Reco data/MC check");

              c.cd(4);
              prepPad();
              hTruthOverData->GetXaxis()->SetRangeUser(10.0, 35.0);
              hTruthOverData->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
              hTruthOverData->GetYaxis()->SetTitle("SIM truth / DATA reco");
              setRatioRange(hTruthOverData, 0.5, 2.5);
              hTruthOverData->Draw("E1");
              {
                TLine l1(10.0, 1.0, 35.0, 1.0);
                l1.SetLineStyle(2);
                l1.SetLineWidth(2);
                l1.Draw("same");
              }
              drawPadTitle("Overall truth/data gap");

              c.cd();
              drawCanvasTitle("Photon master QA: yields and key ratios", 0.028);
              drawCanvasNote(
                {
                  "Displayed bins: 10-12, 12-14, 14-16, 16-18, 18-20, 20-22, 22-24, 24-26, 26-35 GeV",
                  "SIM truth is mapped onto the reco analysis axis when needed"
                },
                0.97, 0.05, 31, 0.020
              );
              SaveCanvas(c, JoinPath(qaMasterDir, "pho_master_yields_and_ratios.png"));
            }

            if (hFakeFracSim && hPuritySim && hMissFracSim && hEffSim && hTruthOverRecoSim && hAfterOverBefore)
            {
              TCanvas c("c_pho_budgetQA", "c_pho_budgetQA", 1700, 1050);
              c.Divide(3, 2, 0.001, 0.001);

              vector<pair<TH1*, string>> pads =
              {
                { hFakeFracSim,     "Fake fraction" },
                { hPuritySim,       "Purity = 1 - fake/reco" },
                { hMissFracSim,     "Miss fraction" },
                { hEffSim,          "Efficiency = 1 - miss/truth" },
                { hTruthOverRecoSim,"SIM truth / SIM reco" },
                { hAfterOverBefore, "DATA after / before unfolding" }
              };

              for (int ipad = 0; ipad < 6; ++ipad)
              {
                c.cd(ipad + 1);
                prepPad();

                TH1* h = pads[(std::size_t)ipad].first;
                if (!h)
                {
                  drawMissingPad("Value", "MISSING");
                  continue;
                }

                h->GetXaxis()->SetRangeUser(10.0, 35.0);
                h->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                h->GetYaxis()->SetTitle("Value");

                if (ipad < 4)
                {
                  h->GetYaxis()->SetRangeUser(0.0, 1.2);
                }
                else
                {
                  setRatioRange(h, 0.5, 2.5);
                }

                h->Draw("E1");

                if (ipad != 0 && ipad != 2)
                {
                  TLine l1(10.0, 1.0, 35.0, 1.0);
                  l1.SetLineStyle(2);
                  l1.SetLineWidth(2);
                  l1.Draw("same");
                }

                drawPadTitle(pads[(std::size_t)ipad].second);
              }

              c.cd();
              drawCanvasTitle("Photon correction-budget QA", 0.028);
              drawCanvasNote(
                {
                  "Top row: reco fake contamination and truth-matching loss terms",
                  "Bottom-right compares the observed data correction to the SIM truth/reco lever-arm"
                },
                0.97, 0.05, 31, 0.020
              );
              SaveCanvas(c, JoinPath(qaBudgetDir, "pho_correction_budget_summary.png"));
            }

            if (hPhoRespSim)
            {
              TH2* hRespTruthNorm = normalizeResponseByTruth(hPhoRespSim, "hPhoResp_truthNorm_QA");
              TH2* hRespRecoNorm  = normalizeResponseByReco (hPhoRespSim, "hPhoResp_recoNorm_QA");
              TH1D* hRespMean     = buildResponseMeanRecoOverTruth(hPhoRespSim, "hPhoResp_meanRecoOverTruth_QA");
              TH1D* hRespWidth    = buildResponseWidthRecoOverTruth(hPhoRespSim, "hPhoResp_widthRecoOverTruth_QA");
              TH1D* hRespDiagFrac = buildResponseDiagonalFraction(hPhoRespSim, "hPhoResp_diagFrac_QA");

              styleHist(hRespMean,     kBlack,    20);
              styleHist(hRespWidth,    kBlue + 1, 20);
              styleHist(hRespDiagFrac, kRed + 1,  20);

              TCanvas c("c_pho_responseShapeQA", "c_pho_responseShapeQA", 1750, 1050);
              c.Divide(3, 2, 0.001, 0.001);

              c.cd(1);
              prepPad();
              if (gPad) gPad->SetLogz(1);
              TH2* hRawDraw = CloneTH2(hPhoRespSim, "hPhoResp_rawDraw_QA");
              if (hRawDraw)
              {
                hRawDraw->SetDirectory(nullptr);
                hRawDraw->SetTitle("");
                hRawDraw->GetXaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
                hRawDraw->GetYaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                double minPos = 1e99;
                for (int ix = 1; ix <= hRawDraw->GetNbinsX(); ++ix)
                {
                  for (int iy = 1; iy <= hRawDraw->GetNbinsY(); ++iy)
                  {
                    const double z = hRawDraw->GetBinContent(ix, iy);
                    if (std::isfinite(z) && z > 0.0 && z < minPos) minPos = z;
                  }
                }
                if (minPos < 1e98) hRawDraw->SetMinimum(std::max(0.5 * minPos, 1e-6));
                hRawDraw->Draw("colz");
                drawPadTitle("Raw response counts");
              }

              c.cd(2);
              prepPad();
              if (hRespTruthNorm)
              {
                hRespTruthNorm->SetTitle("");
                hRespTruthNorm->GetXaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
                hRespTruthNorm->GetYaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                hRespTruthNorm->SetMinimum(0.0);
                hRespTruthNorm->SetMaximum(1.0);
                hRespTruthNorm->Draw("colz");
                drawPadTitle("Truth-bin normalized: P(reco|truth)");
              }
              else
              {
                drawMissingPad("Probability", "MISSING");
              }

              c.cd(3);
              prepPad();
              if (hRespRecoNorm)
              {
                hRespRecoNorm->SetTitle("");
                hRespRecoNorm->GetXaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
                hRespRecoNorm->GetYaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                hRespRecoNorm->SetMinimum(0.0);
                hRespRecoNorm->SetMaximum(1.0);
                hRespRecoNorm->Draw("colz");
                drawPadTitle("Reco-bin normalized: P(truth|reco)");
              }
              else
              {
                drawMissingPad("Probability", "MISSING");
              }

              c.cd(4);
              prepPad();
              if (hRespMean)
              {
                hRespMean->GetXaxis()->SetRangeUser(5.0, 40.0);
                hRespMean->GetYaxis()->SetRangeUser(0.6, 1.4);
                hRespMean->Draw("E1");
                TLine l1(5.0, 1.0, 40.0, 1.0);
                l1.SetLineStyle(2);
                l1.SetLineWidth(2);
                l1.Draw("same");
                drawPadTitle("Mean reco/truth vs truth p_{T}");
              }
              else
              {
                drawMissingPad("Mean(reco / truth)", "MISSING");
              }

              c.cd(5);
              prepPad();
              if (hRespWidth)
              {
                hRespWidth->GetXaxis()->SetRangeUser(5.0, 40.0);
                hRespWidth->GetYaxis()->SetRangeUser(0.0, 0.6);
                hRespWidth->Draw("E1");
                drawPadTitle("RMS(reco/truth) vs truth p_{T}");
              }
              else
              {
                drawMissingPad("RMS(reco / truth)", "MISSING");
              }

              c.cd(6);
              prepPad();
              if (hRespDiagFrac)
              {
                hRespDiagFrac->GetXaxis()->SetRangeUser(5.0, 40.0);
                hRespDiagFrac->GetYaxis()->SetRangeUser(0.0, 1.05);
                hRespDiagFrac->Draw("E1");
                drawPadTitle("Same-bin fraction vs truth p_{T}");
              }
              else
              {
                drawMissingPad("Same-bin fraction", "MISSING");
              }

              c.cd();
              drawCanvasTitle("Photon response-shape QA", 0.028);
              drawCanvasNote(
                {
                  "Top row shows the raw response and its truth/reco normalizations",
                  "Bottom row compresses the matrix into migration bias, width, and diagonal fraction"
                },
                0.97, 0.05, 31, 0.020
              );
              SaveCanvas(c, JoinPath(qaResponseDir, "pho_responseShape_summary.png"));

              if (hRawDraw) delete hRawDraw;
              if (hRespTruthNorm) delete hRespTruthNorm;
              if (hRespRecoNorm) delete hRespRecoNorm;
              if (hRespMean) delete hRespMean;
              if (hRespWidth) delete hRespWidth;
              if (hRespDiagFrac) delete hRespDiagFrac;
            }

            if (hPhoResp_measXtruth)
            {
              TH2* hFeedIn = normalizeRecoFeedIn(hPhoResp_measXtruth, "hPhoResp_recoFeedIn_QA");
              TH1* hFeed58 = CloneTH1(hPhoRecoData, "hPho_feedFromTruth5to8_QA");
              TH1* hFeed810 = CloneTH1(hPhoRecoData, "hPho_feedFromTruth8to10_QA");
              TH1* hFeedSupport = CloneTH1(hPhoRecoData, "hPho_feedFromTruthSupport_QA");

              if (hFeed58)     { hFeed58->SetDirectory(nullptr);     EnsureSumw2(hFeed58);     hFeed58->Reset("ICES"); }
              if (hFeed810)    { hFeed810->SetDirectory(nullptr);    EnsureSumw2(hFeed810);    hFeed810->Reset("ICES"); }
              if (hFeedSupport){ hFeedSupport->SetDirectory(nullptr);EnsureSumw2(hFeedSupport);hFeedSupport->Reset("ICES"); }

              const int iy58  = (hFeedIn ? findExactBinByEdges(hFeedIn->GetYaxis(), 5.0, 8.0)  : -1);
              const int iy810 = (hFeedIn ? findExactBinByEdges(hFeedIn->GetYaxis(), 8.0, 10.0) : -1);

              if (hFeedIn && hFeed58 && hFeed810 && hFeedSupport)
              {
                for (int ix = 1; ix <= hFeedIn->GetNbinsX(); ++ix)
                {
                  const double recoLo = hFeedIn->GetXaxis()->GetBinLowEdge(ix);
                  const double recoHi = hFeedIn->GetXaxis()->GetBinUpEdge(ix);
                  const double recoCtr = hFeedIn->GetXaxis()->GetBinCenter(ix);
                  const int ibRef = hFeed58->GetXaxis()->FindBin(recoCtr);
                  if (ibRef < 1 || ibRef > hFeed58->GetNbinsX()) continue;

                  if (recoLo < 10.0 || recoHi > 35.0) continue;

                  double f58  = (iy58  > 0 ? hFeedIn->GetBinContent(ix, iy58)  : 0.0);
                  double f810 = (iy810 > 0 ? hFeedIn->GetBinContent(ix, iy810) : 0.0);
                  const double fSup = f58 + f810;

                  hFeed58->SetBinContent(ibRef, f58);
                  hFeed810->SetBinContent(ibRef, f810);
                  hFeedSupport->SetBinContent(ibRef, fSup);
                  hFeed58->SetBinError(ibRef, 0.0);
                  hFeed810->SetBinError(ibRef, 0.0);
                  hFeedSupport->SetBinError(ibRef, 0.0);
                }

                styleHist(hFeed58,      kBlue + 1, 24);
                styleHist(hFeed810,     kRed + 1,  25);
                styleHist(hFeedSupport, kBlack,    20);

                TCanvas c("c_pho_supportQA", "c_pho_supportQA", 1600, 760);
                c.Divide(2, 1, 0.001, 0.001);

                c.cd(1);
                prepPad();
                if (hFeedIn)
                {
                  hFeedIn->SetTitle("");
                  hFeedIn->GetXaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                  hFeedIn->GetYaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
                  hFeedIn->GetXaxis()->SetRangeUser(10.0, 35.0);
                  hFeedIn->GetYaxis()->SetRangeUser(5.0, 40.0);
                  hFeedIn->SetMinimum(0.0);
                  hFeedIn->SetMaximum(1.0);
                  hFeedIn->Draw("colz");

                  if (iy58 > 0)
                  {
                    TLine l(10.0, hFeedIn->GetYaxis()->GetBinUpEdge(iy58), 35.0, hFeedIn->GetYaxis()->GetBinUpEdge(iy58));
                    l.SetLineStyle(2);
                    l.SetLineWidth(2);
                    l.Draw("same");
                  }
                  if (iy810 > 0)
                  {
                    TLine l(10.0, hFeedIn->GetYaxis()->GetBinUpEdge(iy810), 35.0, hFeedIn->GetYaxis()->GetBinUpEdge(iy810));
                    l.SetLineStyle(2);
                    l.SetLineWidth(2);
                    l.Draw("same");
                  }

                  drawPadTitle("Reco-bin feed-in from truth bins");
                }
                else
                {
                  drawMissingPad("Fraction", "MISSING");
                }

                c.cd(2);
                prepPad();
                if (hFeed58 && hFeed810 && hFeedSupport)
                {
                  double maxY = 0.0;
                  auto scan = [&](TH1* h)->void
                  {
                    for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                    {
                      const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                      const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                      if (lo < 10.0 || hi > 35.0) continue;
                      const double v = h->GetBinContent(ib) + h->GetBinError(ib);
                      if (std::isfinite(v) && v > maxY) maxY = v;
                    }
                  };
                  scan(hFeed58);
                  scan(hFeed810);
                  scan(hFeedSupport);

                  hFeedSupport->GetXaxis()->SetRangeUser(10.0, 35.0);
                  hFeedSupport->GetXaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                  hFeedSupport->GetYaxis()->SetTitle("Feed-in fraction");
                  hFeedSupport->GetYaxis()->SetRangeUser(0.0, (maxY > 0.0) ? (1.15 * maxY) : 0.20);
                  hFeedSupport->Draw("E1");
                  hFeed58->GetXaxis()->SetRangeUser(10.0, 35.0);
                  hFeed810->GetXaxis()->SetRangeUser(10.0, 35.0);
                  hFeed58->Draw("E1 same");
                  hFeed810->Draw("E1 same");
                  drawPadTitle("Low-p_{T}^{truth} support feed-in");

                  TLegend leg(0.48, 0.68, 0.90, 0.88);
                  leg.SetBorderSize(0);
                  leg.SetFillStyle(0);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.040);
                  leg.AddEntry(hFeedSupport, "Total support feed-in (5-10 GeV truth)", "pe");
                  leg.AddEntry(hFeed58,      "5-8 GeV truth", "pe");
                  leg.AddEntry(hFeed810,     "8-10 GeV truth", "pe");
                  leg.Draw();
                }
                else
                {
                  drawMissingPad("Feed-in fraction", "MISSING");
                }

                c.cd();
                drawCanvasTitle("Photon support-bin feed-in QA", 0.028);
                drawCanvasNote(
                  {
                    "Left: each reco analysis bin normalized to unit truth-bin feed-in",
                    "Right: isolates how much of each displayed reco bin comes from the low-p_{T}^{truth} support bins"
                  },
                  0.97, 0.05, 31, 0.020
                );
                SaveCanvas(c, JoinPath(qaSupportDir, "pho_supportFeedin_summary.png"));
              }

              if (hFeedIn) delete hFeedIn;
              if (hFeed58) delete hFeed58;
              if (hFeed810) delete hFeed810;
              if (hFeedSupport) delete hFeedSupport;
            }

            if (hAfterOverBefore) delete hAfterOverBefore;
            if (hEffSim) delete hEffSim;
            if (hMissFracSim) delete hMissFracSim;
            if (hPuritySim) delete hPuritySim;
            if (hFakeFracSim) delete hFakeFracSim;
            if (hTruthOverData) delete hTruthOverData;
            if (hRecoOverData) delete hRecoOverData;
            if (hTruthOverRecoSim) delete hTruthOverRecoSim;
            if (hUnfoldDataQA) delete hUnfoldDataQA;
            if (hMissesSimQA) delete hMissesSimQA;
            if (hTruthSimQA) delete hTruthSimQA;
            if (hRecoSimQA) delete hRecoSimQA;
            if (hRecoDataQA) delete hRecoDataQA;
          }
        }

        // Photon summary for unfolding RECO pT bins
        vector<string> phoSummary;
        phoSummary.push_back("Photon unfolding summary (PP DATA unfolded to truth pTgamma)");
        phoSummary.push_back(TString::Format("Method: RooUnfoldBayes, iterations=%d, error=kCovToy, Ntoys=%d", kBayesIterPho, kNToysPhoFinal).Data());
        phoSummary.push_back("");

        {
            const auto& analysisRecoBins = UnfoldAnalysisRecoPtBins();
            const int nPtUnf = (int)analysisRecoBins.size();
            for (int i = 0; i < nPtUnf; ++i)
            {
              const PtBin& b = analysisRecoBins[i];
              const double cen = 0.5 * (b.lo + b.hi);

              const int ibTruth = (hPhoUnfoldTruth ? hPhoUnfoldTruth->GetXaxis()->FindBin(cen) : -1);

              const string labCanon = TString::Format("%d-%d GeV", b.lo, b.hi).Data();
              const string labTruth = (hPhoUnfoldTruth ? AxisBinLabel(hPhoUnfoldTruth->GetXaxis(), ibTruth, "GeV", 0) : "N/A");

              const double N = (hPhoUnfoldTruth ? hPhoUnfoldTruth->GetBinContent(ibTruth) : 0.0);
              const double E = (hPhoUnfoldTruth ? hPhoUnfoldTruth->GetBinError  (ibTruth) : 0.0);

              phoSummary.push_back(
                TString::Format("pT^gamma analysis=%s  -> truthBin=%s  N_gamma(unf)=%.6g ± %.6g",
                  labCanon.c_str(), labTruth.c_str(), N, E
                ).Data()
              );
            }
            phoSummary.push_back("");
            phoSummary.push_back("NOTE: only the 9 analysis bins are summarized/plotted: 10-12, 12-14, 14-16, 16-18, 18-20, 20-22, 22-24, 24-26, 26-35 GeV.");
            phoSummary.push_back("      The reco 8-10 GeV bin is retained only as an unfolding support underflow bin.");
            phoSummary.push_back("      The truth 5-8 GeV and 8-10 GeV bins are retained only as truth underflow support bins.");
            phoSummary.push_back("      The reco/truth 35-40 GeV bin is retained only as an overflow support bin.");
            phoSummary.push_back("      RooUnfold still uses the full YAML unfolding axes; only the displayed analysis-bin outputs are restricted here.");
        }

      WriteTextFile(JoinPath(phoDir, "summary_photon_unfolding_bins.txt"), phoSummary);

      // Save photon ROOT outputs
      {
        const string outRoot = JoinPath(phoDir, "rooUnfold_photons.root");
        TFile f(outRoot.c_str(), "RECREATE");
        if (f.IsOpen())
        {
          if (hPhoRespSim)          hPhoRespSim->Write("h2_phoResp_truthVsReco");
          if (hPhoResp_measXtruth)  hPhoResp_measXtruth->Write("h2_phoResp_recoVsTruth");
          if (hPhoRecoData)         hPhoRecoData->Write("h_phoReco_data");
          if (hPhoRecoSim)          hPhoRecoSim->Write("h_phoReco_sim");
          if (hPhoTruthSim)         hPhoTruthSim->Write("h_phoTruth_sim");
          if (hPhoUnfoldTruth)      hPhoUnfoldTruth->Write("h_phoTruth_unfolded_data");
          if (hPhoUnfoldTruth_cov)  hPhoUnfoldTruth_cov->Write("h_phoTruth_unfolded_data_covariance");
          f.Close();
        }
      }

      // ----------------------------------------------------------------------
      // NEW (Step A validation): 1D photon unfolding QA package
      //
      // Output folder:
      //   <phoDir>/validation/
      //
      // Includes:
      //   (1) Iteration stability: rel(stat) + rel(change it vs it-1)
      //   (2) Closure: unfold(reco SIM) -> truth SIM
      //   (3) Half-closure: train(A) unfold(B) using binomial split of response
      // ----------------------------------------------------------------------
      {
        const string phoValDir = JoinPath(phoDir, "validation");
        EnsureDir(phoValDir);

          // -----------------------------
          // (1) Iteration stability (DATA)
          // -----------------------------
          if (hPhoRecoData && hPhoTruthSim && hPhoRecoSim && hPhoResp_measXtruth)
          {
            const int kMaxIt = kMaxBayesIterPhoScan;

            const vector<double>& xIt       = phoIterScan.xIt;
            const vector<double>& exIt      = phoIterScan.exIt;
            const vector<double>& yRelStat  = phoIterScan.yRelStat;
            const vector<double>& eyRelStat = phoIterScan.eyRelStat;
            const vector<double>& yRelDev   = phoIterScan.yRelChange;
            const vector<double>& eyRelDev  = phoIterScan.eyRelChange;
            const vector<double>& yQuad     = phoIterScan.yQuad;
            const vector<double>& eyQuad    = phoIterScan.eyQuad;
            const int bestIt                = phoIterScan.bestIt;
            const double bestQuad           = phoIterScan.bestQuad;

              if ((int)xIt.size() > 0)
              {
                vector<double> xItNo1;
                vector<double> exItNo1;
                vector<double> yRelStatNo1;
                vector<double> eyRelStatNo1;
                vector<double> yRelDevNo1;
                vector<double> eyRelDevNo1;
                vector<double> yQuadNo1;
                vector<double> eyQuadNo1;

                for (size_t i = 0; i < xIt.size(); ++i)
                {
                  if (xIt[i] <= 1.0) continue;
                  xItNo1.push_back(xIt[i]);
                  exItNo1.push_back(exIt[i]);
                  yRelStatNo1.push_back(yRelStat[i]);
                  eyRelStatNo1.push_back(eyRelStat[i]);
                  yRelDevNo1.push_back(yRelDev[i]);
                  eyRelDevNo1.push_back(eyRelDev[i]);
                  yQuadNo1.push_back(yQuad[i]);
                  eyQuadNo1.push_back(eyQuad[i]);
                }

                TCanvas cSt("c_pho_iterStability_relChange_relStat", "c_pho_iterStability_relChange_relStat", 900, 700);
                ApplyCanvasMargins1D(cSt);

                // Determine y-range from content
                double ymax = 0.0;
                for (size_t i = 0; i < yRelStat.size(); ++i) ymax = std::max(ymax, yRelStat[i]);
                for (size_t i = 0; i < yRelDev.size();  ++i) ymax = std::max(ymax, yRelDev[i]);
                if (ymax <= 0.0) ymax = 0.1;
                ymax *= 1.25;

                TH1F frame("frame_phoIt", "", 1, 1.0, (double)kMaxIt + 0.5);
                frame.SetMinimum(0.0);
                frame.SetMaximum(ymax);
                frame.SetTitle("");
                frame.GetXaxis()->SetTitle("Iteration");
                frame.GetYaxis()->SetTitle("Relative quantity");
                frame.Draw("axis");

                // total relative stat. uncertainty (RED)
                TGraphErrors gStat((int)xIt.size(), &xIt[0], &yRelStat[0], &exIt[0], &eyRelStat[0]);
                gStat.SetMarkerStyle(24);
                gStat.SetMarkerSize(1.1);
                gStat.SetMarkerColor(kRed + 1);
                gStat.SetLineColor(kRed + 1);
                gStat.SetLineWidth(2);
                gStat.Draw("P same");

                // total relative deviation (it vs it-1, with it=1 using 0->1 baseline) (BLUE)
                TGraphErrors gDev((int)xIt.size(), &xIt[0], &yRelDev[0], &exIt[0], &eyRelDev[0]);
                gDev.SetMarkerStyle(24);
                gDev.SetMarkerSize(1.1);
                gDev.SetMarkerColor(kBlue + 1);
                gDev.SetLineColor(kBlue + 1);
                gDev.SetLineWidth(2);
                gDev.Draw("P same");

                TLegend leg(0.55, 0.78, 0.89, 0.90);
                leg.SetBorderSize(0);
                leg.SetFillStyle(0);
                leg.SetTextFont(42);
                leg.SetTextSize(0.033);
                leg.AddEntry(&gStat, "total relative stat. uncertainty", "p");
                leg.AddEntry(&gDev,  "total relative deviation (it vs it-1)", "p");
                leg.Draw();

                {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(22);
                    tx.SetTextSize(0.040);
                    tx.DrawLatex(0.50, 0.965, "Iteration Stability, Photon 4 + MBD NS #geq 1, Run24pp");
                }

                {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(13);
                    tx.SetTextSize(0.032);
                    tx.DrawLatex(0.55, 0.74, "1D photon-yield unfolding");
                }

                SaveCanvas(cSt, JoinPath(phoValDir, "pho_unfold_iterStability_relChange_relStat.png"));

                if (!xItNo1.empty())
                {
                  TCanvas cStNo1("c_pho_iterStability_relChange_relStat_noIter1", "c_pho_iterStability_relChange_relStat_noIter1", 900, 700);
                  ApplyCanvasMargins1D(cStNo1);

                  double ymaxNo1 = 0.0;
                  for (size_t i = 0; i < yRelStatNo1.size(); ++i) ymaxNo1 = std::max(ymaxNo1, yRelStatNo1[i]);
                  for (size_t i = 0; i < yRelDevNo1.size();  ++i) ymaxNo1 = std::max(ymaxNo1, yRelDevNo1[i]);
                  if (ymaxNo1 <= 0.0) ymaxNo1 = 0.1;
                  ymaxNo1 *= 1.25;

                  TH1F frameNo1("frame_phoIt_noIter1", "", 1, 1.0, (double)kMaxIt + 0.5);
                  frameNo1.SetMinimum(0.0);
                  frameNo1.SetMaximum(ymaxNo1);
                  frameNo1.SetTitle("");
                  frameNo1.GetXaxis()->SetTitle("Iteration");
                  frameNo1.GetYaxis()->SetTitle("Relative quantity");
                  frameNo1.Draw("axis");

                  TGraphErrors gStatNo1((int)xItNo1.size(), &xItNo1[0], &yRelStatNo1[0], &exItNo1[0], &eyRelStatNo1[0]);
                  gStatNo1.SetMarkerStyle(24);
                  gStatNo1.SetMarkerSize(1.1);
                  gStatNo1.SetMarkerColor(kRed + 1);
                  gStatNo1.SetLineColor(kRed + 1);
                  gStatNo1.SetLineWidth(2);
                  gStatNo1.Draw("P same");

                  TGraphErrors gDevNo1((int)xItNo1.size(), &xItNo1[0], &yRelDevNo1[0], &exItNo1[0], &eyRelDevNo1[0]);
                  gDevNo1.SetMarkerStyle(24);
                  gDevNo1.SetMarkerSize(1.1);
                  gDevNo1.SetMarkerColor(kBlue + 1);
                  gDevNo1.SetLineColor(kBlue + 1);
                  gDevNo1.SetLineWidth(2);
                  gDevNo1.Draw("P same");

                  TLegend legNo1(0.55, 0.78, 0.89, 0.90);
                  legNo1.SetBorderSize(0);
                  legNo1.SetFillStyle(0);
                  legNo1.SetTextFont(42);
                  legNo1.SetTextSize(0.033);
                  legNo1.AddEntry(&gStatNo1, "total relative stat. uncertainty", "p");
                  legNo1.AddEntry(&gDevNo1,  "total relative deviation (it vs it-1)", "p");
                  legNo1.Draw();

                  {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(22);
                      tx.SetTextSize(0.040);
                      tx.DrawLatex(0.50, 0.965, "Iteration Stability, Photon 4 + MBD NS #geq 1, Run24pp");
                  }

                  {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(13);
                      tx.SetTextSize(0.032);
                      tx.DrawLatex(0.55, 0.74, "1D photon-yield unfolding (it #geq 2 only)");
                  }

                  SaveCanvas(cStNo1, JoinPath(phoValDir, "pho_unfold_iterStability_relChange_relStat_noIter1.png"));
                }

                TCanvas cQuad("c_pho_iterStability_quadratureSum", "c_pho_iterStability_quadratureSum", 900, 700);
                ApplyCanvasMargins1D(cQuad);

                double ymaxQ = 0.0;
                for (size_t i = 0; i < yQuad.size(); ++i) ymaxQ = std::max(ymaxQ, yQuad[i]);
                if (ymaxQ <= 0.0) ymaxQ = 0.1;
                ymaxQ *= 1.25;

                TH1F frameQ("frame_phoItQuad", "", 1, 1.0, (double)kMaxIt + 0.5);
                frameQ.SetMinimum(0.0);
                frameQ.SetMaximum(ymaxQ);
                frameQ.SetTitle("");
                frameQ.GetXaxis()->SetTitle("Iteration");
                frameQ.GetYaxis()->SetTitle("Quadrature sum");
                frameQ.Draw("axis");

                TGraphErrors gQuad((int)xIt.size(), &xIt[0], &yQuad[0], &exIt[0], &eyQuad[0]);
                gQuad.SetMarkerStyle(20);
                gQuad.SetMarkerSize(1.1);
                gQuad.SetMarkerColor(kBlack);
                gQuad.SetLineColor(kBlack);
                gQuad.SetLineWidth(2);
                gQuad.Draw("LP same");

                TGraphErrors gBest;
                if (bestIt > 0 && std::isfinite(bestQuad))
                {
                  const double bestX[1]  = { (double)bestIt };
                  const double bestY[1]  = { bestQuad };
                  const double bestEX[1] = { 0.0 };
                  const double bestEY[1] = { 0.0 };

                  gBest = TGraphErrors(1, bestX, bestY, bestEX, bestEY);
                  gBest.SetMarkerStyle(29);
                  gBest.SetMarkerSize(1.6);
                  gBest.SetMarkerColor(kRed + 1);
                  gBest.SetLineColor(kRed + 1);
                  gBest.Draw("P same");
                }

                TLegend legQ(0.55, 0.78, 0.89, 0.90);
                legQ.SetBorderSize(0);
                legQ.SetFillStyle(0);
                legQ.SetTextFont(42);
                legQ.SetTextSize(0.033);
                legQ.AddEntry(&gQuad, "quadrature sum", "lp");
                if (bestIt > 0 && std::isfinite(bestQuad))
                {
                  legQ.AddEntry(&gBest, TString::Format("minimum: it=%d", bestIt).Data(), "p");
                }
                legQ.Draw();

                {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(22);
                    tx.SetTextSize(0.040);
                    tx.DrawLatex(0.50, 0.965, "Quadrature-sum optimization, Photon 4 + MBD NS #geq 1, Run24pp");
                }

                {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(13);
                    tx.SetTextSize(0.032);
                    tx.DrawLatex(0.55, 0.74, "1D photon-yield unfolding");
                }

                SaveCanvas(cQuad, JoinPath(phoValDir, "pho_unfold_iterStability_quadratureSum.png"));

                if (!xItNo1.empty())
                {
                  TCanvas cQuadNo1("c_pho_iterStability_quadratureSum_noIter1", "c_pho_iterStability_quadratureSum_noIter1", 900, 700);
                  ApplyCanvasMargins1D(cQuadNo1);

                  double ymaxQNo1 = 0.0;
                  for (size_t i = 0; i < yQuadNo1.size(); ++i) ymaxQNo1 = std::max(ymaxQNo1, yQuadNo1[i]);
                  if (ymaxQNo1 <= 0.0) ymaxQNo1 = 0.1;
                  ymaxQNo1 *= 1.25;

                  TH1F frameQNo1("frame_phoItQuad_noIter1", "", 1, 1.0, (double)kMaxIt + 0.5);
                  frameQNo1.SetMinimum(0.0);
                  frameQNo1.SetMaximum(ymaxQNo1);
                  frameQNo1.SetTitle("");
                  frameQNo1.GetXaxis()->SetTitle("Iteration");
                  frameQNo1.GetYaxis()->SetTitle("Quadrature sum");
                  frameQNo1.Draw("axis");

                  TGraphErrors gQuadNo1((int)xItNo1.size(), &xItNo1[0], &yQuadNo1[0], &exItNo1[0], &eyQuadNo1[0]);
                  gQuadNo1.SetMarkerStyle(20);
                  gQuadNo1.SetMarkerSize(1.1);
                  gQuadNo1.SetMarkerColor(kBlack);
                  gQuadNo1.SetLineColor(kBlack);
                  gQuadNo1.SetLineWidth(2);
                  gQuadNo1.Draw("LP same");

                  TLegend legQNo1(0.55, 0.78, 0.89, 0.90);
                  legQNo1.SetBorderSize(0);
                  legQNo1.SetFillStyle(0);
                  legQNo1.SetTextFont(42);
                  legQNo1.SetTextSize(0.033);
                  legQNo1.AddEntry(&gQuadNo1, "quadrature sum", "lp");
                  legQNo1.Draw();

                  {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(22);
                      tx.SetTextSize(0.040);
                      tx.DrawLatex(0.50, 0.965, "Quadrature-sum optimization, Photon 4 + MBD NS #geq 1, Run24pp");
                  }

                  {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(13);
                      tx.SetTextSize(0.032);
                      tx.DrawLatex(0.55, 0.74, "1D photon-yield unfolding (it #geq 2 only)");
                  }

                  SaveCanvas(cQuadNo1, JoinPath(phoValDir, "pho_unfold_iterStability_quadratureSum_noIter1.png"));
                }

                vector<string> iterSummary;
                iterSummary.push_back("Photon iteration-stability summary");
                iterSummary.push_back(TString::Format("best iteration from quadrature sum = %d", bestIt).Data());
                iterSummary.push_back(TString::Format("minimum quadrature sum = %.10g", bestQuad).Data());
                iterSummary.push_back("relChange is defined as the successive-iterate difference it vs (it-1), with iteration 1 using the explicit 0->1 baseline.");
                iterSummary.push_back("iteration 0 baseline is the measured input histogram mapped onto the unfolding comparison axis, so iteration 1 is included in the best-iteration choice.");
                iterSummary.push_back("additional plots with suffix _noIter1 exclude iteration 1 from the display only; they do not change the scan or best-iteration choice.");
                iterSummary.push_back("relStat is built from the full RooUnfold::kCovToy covariance restricted to the 9 analysis p_{T}^{#gamma} bins.");
                iterSummary.push_back("support / underflow / overflow bins are excluded from the stability scalar.");
                iterSummary.push_back("quadrature sum definition: sqrt(relStat^2 + relChange^2)");
                iterSummary.push_back("");

                for (size_t i = 0; i < xIt.size(); ++i)
                {
                  iterSummary.push_back(
                    TString::Format("it=%d  relStat=%.10g  relChange=%.10g  quadratureSum=%.10g",
                                    (int)xIt[i], yRelStat[i], yRelDev[i], yQuad[i]).Data()
                  );
                }

                WriteTextFile(JoinPath(phoValDir, "pho_unfold_iterStability_bestIteration.txt"), iterSummary);

                cout << ANSI_BOLD_CYN
                     << "[PHO ITER QA] best iteration from quadrature sum = " << bestIt
                     << "  minimum = " << bestQuad
                     << ANSI_RESET << "\n";
              }
          }

          // -----------------------------
          // (2) Closure: unfold(reco SIM) -> truth SIM
          // -----------------------------
          if (hPhoResp_measXtruth)
          {
            TH1D* hPhoRecoClosure = hPhoResp_measXtruth->ProjectionX(
              "h_phoRecoClosure_fromResponse",
              0, hPhoResp_measXtruth->GetNbinsY() + 1, "e"
            );
            TH1D* hPhoTruthClosure = hPhoResp_measXtruth->ProjectionY(
              "h_phoTruthClosure_fromResponse",
              0, hPhoResp_measXtruth->GetNbinsX() + 1, "e"
            );

            if (hPhoRecoClosure)  { hPhoRecoClosure->SetDirectory(nullptr);  EnsureSumw2(hPhoRecoClosure); }
            if (hPhoTruthClosure) { hPhoTruthClosure->SetDirectory(nullptr); EnsureSumw2(hPhoTruthClosure); }

            if (hPhoRecoClosure && hPhoTruthClosure)
            {
              RooUnfoldResponse respPhoClosure(
                hPhoRecoClosure, hPhoTruthClosure, hPhoResp_measXtruth,
                "respPhoClosure", "respPhoClosure"
              );

              RooUnfoldBayes uC(&respPhoClosure, hPhoRecoClosure, kBayesIterPho);
              uC.SetVerbose(0);
              uC.SetNToys(kNToysPhoScan);

              TH1* hUnfC = nullptr;
              if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
              hUnfC = uC.Hreco(RooUnfold::kCovToy);
              if (gSystem) gSystem->RedirectOutput(0);
              if (hUnfC)
              {
                hUnfC->SetDirectory(nullptr);
                EnsureSumw2(hUnfC);

                TH1* hRat = CloneTH1(hUnfC, "h_pho_closure_unfoldedOverTruth");
                if (hRat)
                {
                  hRat->SetDirectory(nullptr);
                  EnsureSumw2(hRat);
                  hRat->Divide(hPhoTruthClosure);

                  hRat->SetTitle("");
                  hRat->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                  hRat->GetXaxis()->SetTitleOffset(1.15);
                  hRat->GetYaxis()->SetTitle("Closure: unfolded MC / truth MC");
                  hRat->SetMarkerStyle(24);
                  hRat->SetMarkerSize(1.1);
                  hRat->SetLineWidth(2);

                  TCanvas c("c_pho_closure", "c_pho_closure", 900, 700);
                  ApplyCanvasMargins1D(c);

                  const double xPlotMin = 10.0;
                  const double xPlotMax = 35.0;

                  const double ymin = 0.95;
                  const double ymax = 1.05;

                  hRat->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                  hRat->GetYaxis()->SetRangeUser(ymin, ymax);
                  hRat->Draw("E1");

                  TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                  l1.SetLineStyle(2);
                  l1.SetLineWidth(2);
                  l1.Draw("same");

                    // Top-left Bayes + Step label (match 2D closure styling)
                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(13);
                      tx.SetTextSize(0.038);
                      tx.DrawLatex(0.20, 0.90, TString::Format("Bayes it=%d (photon)", kBayesIterPho).Data());
                      tx.DrawLatex(0.20, 0.855, "1D photon-yield unfolding");
                     }

                     // Centered title
                     {
                        TLatex tx;
                        tx.SetNDC(true);
                        tx.SetTextFont(42);
                        tx.SetTextAlign(22);
                        tx.SetTextSize(0.040);
                        tx.DrawLatex(0.50, 0.965, "Photon closure test, unfold(reco) #rightarrow truth (SIM)");
                    }

                    SaveCanvas(c, JoinPath(phoValDir, "pho_closure_unfoldedOverTruth_vs_pTgamma.png"));

                  delete hRat;
                }

                delete hUnfC;
              }
            }

            if (hPhoRecoClosure) delete hPhoRecoClosure;
            if (hPhoTruthClosure) delete hPhoTruthClosure;
        }

          // -----------------------------
          // (3) Half-closure: train(A) unfold(B)
          //     Split response (reco x truth) bin-by-bin with Binomial(0.5)
          // -----------------------------
          if (hPhoResp_measXtruth && hPhoRecoSim && hPhoTruthSim)
          {
            const int nReco  = hPhoResp_measXtruth->GetNbinsX();
            const int nTruth = hPhoResp_measXtruth->GetNbinsY();

            unsigned int seed = 1337u;
            seed = seed * 131u + (unsigned int)nReco;
            seed = seed * 131u + (unsigned int)nTruth;
            std::mt19937 rng(seed);

            TH2* hRspA = CloneTH2(hPhoResp_measXtruth, "h2_phoResp_recoVsTruth_halfA");
            TH2* hRspB = CloneTH2(hPhoResp_measXtruth, "h2_phoResp_recoVsTruth_halfB");

            TH1* hMeasA  = CloneTH1(hPhoRecoSim,  "h_phoRecoSim_halfA");
            TH1* hMeasB  = CloneTH1(hPhoRecoSim,  "h_phoRecoSim_halfB");
            TH1* hTruthA = CloneTH1(hPhoTruthSim, "h_phoTruthSim_halfA");
            TH1* hTruthB = CloneTH1(hPhoTruthSim, "h_phoTruthSim_halfB");

            if (hRspA)  { hRspA->SetDirectory(nullptr);  hRspA->Reset("ICES");  hRspA->Sumw2(); }
            if (hRspB)  { hRspB->SetDirectory(nullptr);  hRspB->Reset("ICES");  hRspB->Sumw2(); }
            if (hMeasA) { hMeasA->SetDirectory(nullptr); hMeasA->Reset("ICES"); EnsureSumw2(hMeasA); }
            if (hMeasB) { hMeasB->SetDirectory(nullptr); hMeasB->Reset("ICES"); EnsureSumw2(hMeasB); }
            if (hTruthA){ hTruthA->SetDirectory(nullptr);hTruthA->Reset("ICES");EnsureSumw2(hTruthA); }
            if (hTruthB){ hTruthB->SetDirectory(nullptr);hTruthB->Reset("ICES");EnsureSumw2(hTruthB); }

            if (!hRspA || !hRspB || !hMeasA || !hMeasB || !hTruthA || !hTruthB)
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Photon half-closure: failed to allocate split response/marginals. Skipping.\n"
                   << ANSI_RESET;
            }
            else
            {
              vector<double> measA((std::size_t)nReco + 2, 0.0), measB((std::size_t)nReco + 2, 0.0);
              vector<double> truA ((std::size_t)nTruth + 2, 0.0), truB ((std::size_t)nTruth + 2, 0.0);

              for (int ix = 0; ix <= nReco + 1; ++ix)
              {
                for (int iy = 0; iy <= nTruth + 1; ++iy)
                {
                  const double nRaw = hPhoResp_measXtruth->GetBinContent(ix, iy);
                  if (!(nRaw > 0.0))
                  {
                    hRspA->SetBinContent(ix, iy, 0.0); hRspA->SetBinError(ix, iy, 0.0);
                    hRspB->SetBinContent(ix, iy, 0.0); hRspB->SetBinError(ix, iy, 0.0);
                    continue;
                  }

                  const long long N = (long long)std::llround(nRaw);
                  if (N <= 0LL)
                  {
                    hRspA->SetBinContent(ix, iy, 0.0); hRspA->SetBinError(ix, iy, 0.0);
                    hRspB->SetBinContent(ix, iy, 0.0); hRspB->SetBinError(ix, iy, 0.0);
                    continue;
                  }

                  std::binomial_distribution<long long> d(N, 0.5);
                  const long long NA = d(rng);
                  const long long NB = N - NA;

                  hRspA->SetBinContent(ix, iy, (double)NA);
                  hRspA->SetBinError  (ix, iy, std::sqrt((double)NA));

                  hRspB->SetBinContent(ix, iy, (double)NB);
                  hRspB->SetBinError  (ix, iy, std::sqrt((double)NB));

                  if (ix >= 0 && ix <= nReco + 1)
                  {
                    measA[(std::size_t)ix] += (double)NA;
                    measB[(std::size_t)ix] += (double)NB;
                  }
                  if (iy >= 0 && iy <= nTruth + 1)
                  {
                    truA[(std::size_t)iy] += (double)NA;
                    truB[(std::size_t)iy] += (double)NB;
                  }
                }
              }

              for (int ib = 0; ib <= nReco + 1; ++ib)
              {
                hMeasA->SetBinContent(ib, measA[(std::size_t)ib]);
                hMeasA->SetBinError  (ib, std::sqrt(measA[(std::size_t)ib]));
                hMeasB->SetBinContent(ib, measB[(std::size_t)ib]);
                hMeasB->SetBinError  (ib, std::sqrt(measB[(std::size_t)ib]));
              }

              for (int ib = 0; ib <= nTruth + 1; ++ib)
              {
                hTruthA->SetBinContent(ib, truA[(std::size_t)ib]);
                hTruthA->SetBinError  (ib, std::sqrt(truA[(std::size_t)ib]));
                hTruthB->SetBinContent(ib, truB[(std::size_t)ib]);
                hTruthB->SetBinError  (ib, std::sqrt(truB[(std::size_t)ib]));
              }

              RooUnfoldResponse respPhoA(hMeasA, hTruthA, hRspA, "respPho_halfA", "respPho_halfA");

              RooUnfoldBayes uH(&respPhoA, hMeasB, kBayesIterPho);
              uH.SetVerbose(0);
              uH.SetNToys(kNToysPhoScan);

              TH1* hUnfB = nullptr;
              if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
              hUnfB = uH.Hreco(RooUnfold::kCovToy);
              if (gSystem) gSystem->RedirectOutput(0);
              if (hUnfB)
              {
                hUnfB->SetDirectory(nullptr);
                EnsureSumw2(hUnfB);

                TH1* hRat = CloneTH1(hUnfB, "h_pho_halfClosure_unfoldedOverTruth");
                if (hRat)
                {
                  hRat->SetDirectory(nullptr);
                  EnsureSumw2(hRat);
                  hRat->Divide(hTruthB);

                  hRat->SetTitle("");
                  hRat->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                  hRat->GetXaxis()->SetTitleOffset(1.15);
                  hRat->GetYaxis()->SetTitle("Half-closure: unfolded MC / truth MC");
                  hRat->SetMarkerStyle(24);
                  hRat->SetMarkerSize(1.1);
                  hRat->SetLineWidth(2);

                  TCanvas c("c_pho_halfClosure", "c_pho_halfClosure", 900, 700);
                  ApplyCanvasMargins1D(c);

                  const double xPlotMin = 10.0;
                  const double xPlotMax = 35.0;

                  double ymin =  1e99;
                  double ymax = -1e99;
                  bool havePoint = false;

                  for (int ib = 1; ib <= hRat->GetNbinsX(); ++ib)
                  {
                    const double lo = hRat->GetXaxis()->GetBinLowEdge(ib);
                    const double hi = hRat->GetXaxis()->GetBinUpEdge(ib);
                    if (lo < xPlotMin || hi > xPlotMax) continue;

                    const double y  = hRat->GetBinContent(ib);
                    const double ey = hRat->GetBinError(ib);
                    if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                    if (y <= 0.0 && ey <= 0.0) continue;

                    havePoint = true;
                    ymin = std::min(ymin, y - ey);
                    ymax = std::max(ymax, y + ey);
                  }

                  if (!havePoint || !(ymin < ymax))
                  {
                    ymin = 0.90;
                    ymax = 1.10;
                  }
                  else
                  {
                    const double pad = std::max(0.15 * (ymax - ymin), 0.01);
                    ymin -= pad;
                    ymax += pad;

                    if (ymin > 1.0) ymin = 1.0 - pad;
                    if (ymax < 1.0) ymax = 1.0 + pad;
                  }

                  hRat->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                  hRat->GetYaxis()->SetRangeUser(ymin, ymax);
                  hRat->Draw("E1");

                  TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                  l1.SetLineStyle(2);
                  l1.SetLineWidth(2);
                  l1.Draw("same");

                  // Top-left Bayes + Step label (match 2D half-closure styling)
                  {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(13);
                      tx.SetTextSize(0.038);
                      tx.DrawLatex(0.20, 0.90, TString::Format("Bayes it=%d (photon)", kBayesIterPho).Data());
                      tx.DrawLatex(0.20, 0.855, "1D photon-yield unfolding");
                  }

                  // Centered title
                  {
                      TLatex tx;
                      tx.SetNDC(true);
                      tx.SetTextFont(42);
                      tx.SetTextAlign(22);
                      tx.SetTextSize(0.040);
                      tx.DrawLatex(0.50, 0.965, "Photon half-closure test, train(A) unfold(B) (SIM)");
                  }

                  SaveCanvas(c, JoinPath(phoValDir, "pho_halfClosure_unfoldedOverTruth_vs_pTgamma.png"));

                  delete hRat;
                }

                delete hUnfB;
              }
            }

            if (hTruthB) delete hTruthB;
            if (hTruthA) delete hTruthA;
            if (hMeasB)  delete hMeasB;
            if (hMeasA)  delete hMeasA;
            if (hRspB)   delete hRspB;
            if (hRspA)   delete hRspA;
            }
          }

          {
            const string srcLog = JoinPath(phoDir, "pho_unfolded_truth_pTgamma_overlay_logy.png");
            const string dstLog = JoinPath(phoYieldDir, "pho_unfolded_truth_pTgamma_overlay_logy.png");
            const string srcRawPurUnf = JoinPath(phoDir, "pho_raw_purityCorr_unfolded_pTgamma_overlay_logy.png");
            const string dstRawPurUnf = JoinPath(phoYieldDir, "pho_raw_purityCorr_unfolded_pTgamma_overlay_logy.png");

            if (gSystem)
            {
              gSystem->Unlink(dstLog.c_str());
              gSystem->Rename(srcLog.c_str(), dstLog.c_str());

              gSystem->Unlink(dstRawPurUnf.c_str());
              gSystem->Rename(srcRawPurUnf.c_str(), dstRawPurUnf.c_str());
            }
        }

        {
            const string phoDir = JoinPath(outBase, "photons_INCLUSIVE_QA");
            EnsureDir(phoDir);
            const string phoYieldDir = JoinPath(phoDir, "yields");
            EnsureDir(phoYieldDir);

            do
            {
              const string phoRecoName   = "h_unfoldRecoPho_pTgamma_ppg12obj";
              const string phoTruthName  = "h_unfoldTruthPho_pTgamma_ppg12obj";
              const string phoRespName   = "h2_unfoldResponsePho_pTgamma_ppg12obj";
              const string phoFakesName  = "h_unfoldRecoPhoFakes_pTgamma_ppg12obj";       // SIM-only bookkeeping
              const string phoMissesName = "h_unfoldTruthPhoMisses_pTgamma_ppg12obj";     // SIM-only bookkeeping
      
              TH1* hPhoRecoData_in        = GetObj<TH1>(dsData, phoRecoName,   true, true, true);
              TH1* hPhoRecoSim_in         = GetObj<TH1>(dsSim,  phoRecoName,   true, true, true);
              TH1* hPhoTruthSim_in        = GetObj<TH1>(dsSim,  phoTruthName,  true, true, true);
              TH2* hPhoRespSim_in         = GetObj<TH2>(dsSim,  phoRespName,   true, true, true);
              TH1* hPhoRecoFakesSim_in    = GetObj<TH1>(dsSim,  phoFakesName,  true, true, true);
              TH1* hPhoTruthMissesSim_in  = GetObj<TH1>(dsSim,  phoMissesName, true, true, true);
      
              auto MissingWhy = [&](Dataset& ds, const string& relName)->string
              {
                const string fp = FullPath(ds, relName);
                auto it = ds.missingReason.find(fp);
                if (it == ds.missingReason.end()) return "";
                return it->second;
              };
      
              auto PrintGet = [&](Dataset& ds, const string& relName, TObject* obj)
              {
                const string fp = FullPath(ds, relName);
                cout << "  [GET] " << ds.label << " :: " << fp << " : ";
      
                if (!obj)
                {
                  const string why = MissingWhy(ds, relName);
                  cout << ANSI_BOLD_RED << "MISSING" << ANSI_RESET;
                  if (!why.empty()) cout << " (" << why << ")";
                  cout << "\n";
                  return;
                }
      
                cout << ANSI_BOLD_GRN << "FOUND" << ANSI_RESET << " (" << obj->ClassName() << ")";
                if (auto* h2 = dynamic_cast<TH2*>(obj))
                {
                  cout << "  entries=" << h2->GetEntries()
                       << "  nbinsX=" << h2->GetNbinsX()
                       << "  nbinsY=" << h2->GetNbinsY();
                }
                else if (auto* h1 = dynamic_cast<TH1*>(obj))
                {
                  cout << "  entries=" << h1->GetEntries()
                       << "  nbinsX=" << h1->GetNbinsX();
                }
                cout << "\n";
              };
      
              auto DumpTH1Summary = [&](const string& tag, TH1* h)->void
              {
                cout << ANSI_BOLD_CYN << "[DEBUG TH1] " << tag << ANSI_RESET << "\n";
                if (!h)
                {
                  cout << "  <null>\n";
                  return;
                }
      
                cout << "  name=" << h->GetName()
                     << "  title=" << h->GetTitle()
                     << "  nbins=" << h->GetNbinsX()
                     << "  entries=" << h->GetEntries()
                     << "  integral(width)=" << h->Integral(1, h->GetNbinsX(), "width")
                     << "  integral=" << h->Integral(1, h->GetNbinsX())
                     << "\n";
      
                for (int ib = 0; ib <= h->GetNbinsX() + 1; ++ib)
                {
                  const double lo  = (ib >= 1 && ib <= h->GetNbinsX()) ? h->GetXaxis()->GetBinLowEdge(ib) : -999.0;
                  const double hi  = (ib >= 1 && ib <= h->GetNbinsX()) ? h->GetXaxis()->GetBinUpEdge(ib)  : -999.0;
                  const double val = h->GetBinContent(ib);
                  const double err = h->GetBinError(ib);
      
                  if (ib == 0)
                  {
                    cout << "    UF"
                         << "  content=" << val
                         << "  error="   << err << "\n";
                  }
                  else if (ib == h->GetNbinsX() + 1)
                  {
                    cout << "    OF"
                         << "  content=" << val
                         << "  error="   << err << "\n";
                  }
                  else
                  {
                    cout << "    bin " << ib
                         << "  [" << lo << "," << hi << "]"
                         << "  content=" << val
                         << "  error="   << err << "\n";
                  }
                }
              };
      
              auto DumpTH2Summary = [&](const string& tag, TH2* h)->void
              {
                cout << ANSI_BOLD_CYN << "[DEBUG TH2] " << tag << ANSI_RESET << "\n";
                if (!h)
                {
                  cout << "  <null>\n";
                  return;
                }
      
                cout << "  name=" << h->GetName()
                     << "  title=" << h->GetTitle()
                     << "  nbinsX=" << h->GetNbinsX()
                     << "  nbinsY=" << h->GetNbinsY()
                     << "  entries=" << h->GetEntries()
                     << "  integral=" << h->Integral(1, h->GetNbinsX(), 1, h->GetNbinsY())
                     << "\n";
      
                for (int ix = 1; ix <= h->GetNbinsX(); ++ix)
                {
                  double rowInt = 0.0;
                  double rowIntWidth = 0.0;
                  for (int iy = 1; iy <= h->GetNbinsY(); ++iy)
                  {
                    rowInt += h->GetBinContent(ix, iy);
                    rowIntWidth += h->GetBinContent(ix, iy) * h->GetYaxis()->GetBinWidth(iy);
                  }
      
                  cout << "    x-bin " << ix
                       << "  [" << h->GetXaxis()->GetBinLowEdge(ix)
                       << ","   << h->GetXaxis()->GetBinUpEdge(ix) << "]"
                       << "  rowIntegral=" << rowInt
                       << "  rowIntegral(widthY)=" << rowIntWidth
                       << "\n";
                }
              };
      
              auto DumpProjectionYSummary = [&](const string& tag, TH2* h2, int ix)->void
              {
                cout << ANSI_BOLD_CYN << "[DEBUG PROJY] " << tag << ANSI_RESET << "\n";
                if (!h2)
                {
                  cout << "  source TH2 is null\n";
                  return;
                }
                if (ix < 1 || ix > h2->GetNbinsX())
                {
                  cout << "  requested ix=" << ix << " is out of range 1.." << h2->GetNbinsX() << "\n";
                  return;
                }
      
                TH1D* hp = h2->ProjectionY(
                  TString::Format("hDebugProjY_%s_ix%d", h2->GetName(), ix).Data(),
                  ix, ix, "e"
                );
                if (!hp)
                {
                  cout << "  ProjectionY returned null\n";
                  return;
                }
      
                hp->SetDirectory(nullptr);
                EnsureSumw2(hp);
      
                cout << "  x-bin " << ix
                     << "  [" << h2->GetXaxis()->GetBinLowEdge(ix)
                     << ","   << h2->GetXaxis()->GetBinUpEdge(ix) << "]"
                     << "  projIntegral=" << hp->Integral(1, hp->GetNbinsX())
                     << "  projIntegral(width)=" << hp->Integral(1, hp->GetNbinsX(), "width")
                     << "\n";
      
                for (int iy = 0; iy <= hp->GetNbinsX() + 1; ++iy)
                {
                  if (iy == 0)
                  {
                    cout << "    UF"
                         << "  content=" << hp->GetBinContent(iy)
                         << "  error="   << hp->GetBinError(iy) << "\n";
                  }
                  else if (iy == hp->GetNbinsX() + 1)
                  {
                    cout << "    OF"
                         << "  content=" << hp->GetBinContent(iy)
                         << "  error="   << hp->GetBinError(iy) << "\n";
                  }
                  else
                  {
                    cout << "    y-bin " << iy
                         << "  [" << hp->GetXaxis()->GetBinLowEdge(iy)
                         << ","   << hp->GetXaxis()->GetBinUpEdge(iy) << "]"
                         << "  content=" << hp->GetBinContent(iy)
                         << "  error="   << hp->GetBinError(iy) << "\n";
                  }
                }
      
                delete hp;
              };
      
              struct IterScanSummary
              {
                  vector<double> xIt;
                  vector<double> exIt;
                  vector<double> yRelStat;
                  vector<double> eyRelStat;
                  vector<double> yRelChange;
                  vector<double> eyRelChange;
                  vector<double> yQuad;
                  vector<double> eyQuad;
                  int bestIt = -1;
                  double bestQuad = std::numeric_limits<double>::max();
              };
      
                auto BuildPhotonIterScan = [&](RooUnfoldResponse& resp,
                                                 TH1* hRecoMeasured,
                                                 TH1* hTruthTemplate,
                                                 int kMaxIt,
                                                 int nToys)->IterScanSummary
                {
                    IterScanSummary s;
                    if (!hRecoMeasured || !hTruthTemplate) return s;

                    vector<int> selTruthBins;
                    {
                      const auto& analysisBins = UnfoldAnalysisRecoPtBins();
                      TAxis* axTruth = hTruthTemplate->GetXaxis();

                      for (const auto& b : analysisBins)
                      {
                        int ibMatch = -1;
                        for (int ib = 1; ib <= axTruth->GetNbins(); ++ib)
                        {
                          const double lo = axTruth->GetBinLowEdge(ib);
                          const double hi = axTruth->GetBinUpEdge(ib);
                          if (std::fabs(lo - (double)b.lo) < 1e-9 &&
                              std::fabs(hi - (double)b.hi) < 1e-9)
                          {
                            ibMatch = ib;
                            break;
                          }
                        }

                        if (ibMatch > 0) selTruthBins.push_back(ibMatch);
                      }
                    }

                    if (selTruthBins.empty()) return s;

                    TH1* hPrev = CloneTH1(hTruthTemplate, "hPhoIter0Baseline");
                    if (!hPrev) return s;

                    hPrev->SetDirectory(nullptr);
                    EnsureSumw2(hPrev);
                    hPrev->Reset("ICES");

                    for (int ib = 1; ib <= hPrev->GetNbinsX(); ++ib)
                    {
                      const double x = hPrev->GetXaxis()->GetBinCenter(ib);
                      const int iReco = hRecoMeasured->GetXaxis()->FindBin(x);
                      if (iReco < 1 || iReco > hRecoMeasured->GetNbinsX()) continue;

                      hPrev->SetBinContent(ib, hRecoMeasured->GetBinContent(iReco));
                      hPrev->SetBinError  (ib, hRecoMeasured->GetBinError  (iReco));
                    }

                    for (int it = 1; it <= kMaxIt; ++it)
                    {
                      RooUnfoldBayes uIt(&resp, hRecoMeasured, it);
                      uIt.SetVerbose(0);
                      uIt.SetNToys(nToys);

                      TH1* hIt = nullptr;
                      if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
                      hIt = uIt.Hreco(RooUnfold::kCovToy);
                      const auto cov = uIt.Ereco(RooUnfold::kCovToy);
                      if (gSystem) gSystem->RedirectOutput(0);
                      if (!hIt) continue;

                      hIt->SetDirectory(nullptr);
                      EnsureSumw2(hIt);

                      double sumV = 0.0;
                      for (const int ib : selTruthBins)
                      {
                        if (ib < 1 || ib > hIt->GetNbinsX()) continue;
                        sumV += hIt->GetBinContent(ib);
                      }

                      double sumCov = 0.0;
                      const int nCovRows = cov.GetNrows();
                      const int nCovCols = cov.GetNcols();

                      for (size_t i = 0; i < selTruthBins.size(); ++i)
                      {
                        const int ii = selTruthBins[i] - 1;
                        if (ii < 0 || ii >= nCovRows) continue;

                        for (size_t j = 0; j < selTruthBins.size(); ++j)
                        {
                          const int jj = selTruthBins[j] - 1;
                          if (jj < 0 || jj >= nCovCols) continue;
                          sumCov += cov(ii, jj);
                        }
                      }

                      double relStat = 0.0;
                      if (sumV > 0.0) relStat = std::sqrt(std::max(0.0, sumCov)) / sumV;

                      double relDev = 0.0;
                      {
                        double num2 = 0.0;
                        double den2 = 0.0;

                        for (const int ib : selTruthBins)
                        {
                          if (ib < 1 || ib > hIt->GetNbinsX()) continue;

                          const double v  = hIt->GetBinContent(ib);
                          const double vp = hPrev->GetBinContent(ib);
                          const double d  = v - vp;

                          num2 += d * d;
                          den2 += v * v;
                        }

                        if (den2 > 0.0) relDev = std::sqrt(num2 / den2);
                      }

                      const double quad = std::sqrt(relStat * relStat + relDev * relDev);

                      s.xIt.push_back((double)it);
                      s.exIt.push_back(0.0);
                      s.yRelStat.push_back(relStat);
                      s.eyRelStat.push_back(0.0);
                      s.yRelChange.push_back(relDev);
                      s.eyRelChange.push_back(0.0);
                      s.yQuad.push_back(quad);
                      s.eyQuad.push_back(0.0);

                      if (quad < s.bestQuad)
                      {
                        s.bestQuad = quad;
                        s.bestIt = it;
                      }

                      delete hPrev;
                      hPrev = hIt;
                    }

                    if (hPrev) delete hPrev;

                    return s;
                  };
        
                auto BuildXJIterScan = [&](RooUnfoldResponse& resp,
                                           TH1* hMeasuredGlob,
                                           TH1* hTruthTemplateGlob,
                                           TH2* hMeasured2D,
                                           TH2* hTruthTemplate2D,
                                           int kMaxIt,
                                           int nToys,
                                           const string& nameSuffix)->IterScanSummary
                {
                  IterScanSummary s;
                  if (!hMeasuredGlob || !hTruthTemplateGlob || !hMeasured2D || !hTruthTemplate2D) return s;

                  const int nGlobMeasured2D = (hMeasured2D->GetNbinsX() + 2) * (hMeasured2D->GetNbinsY() + 2);
                  const int nGlobTruth2D    = (hTruthTemplate2D->GetNbinsX() + 2) * (hTruthTemplate2D->GetNbinsY() + 2);

                  if (hMeasuredGlob->GetNbinsX() != nGlobMeasured2D) return s;
                  if (hTruthTemplateGlob->GetNbinsX() != nGlobTruth2D) return s;

                  vector<int> selTruthXBins;
                  vector<int> selTruthGlobalBins;

                  {
                    const auto& analysisBins = UnfoldAnalysisRecoPtBins();
                    TAxis* axTruth = hTruthTemplate2D->GetXaxis();

                    for (const auto& b : analysisBins)
                    {
                      int ixMatch = -1;
                      for (int ix = 1; ix <= axTruth->GetNbins(); ++ix)
                      {
                        const double lo = axTruth->GetBinLowEdge(ix);
                        const double hi = axTruth->GetBinUpEdge(ix);
                        if (std::fabs(lo - (double)b.lo) < 1e-9 &&
                            std::fabs(hi - (double)b.hi) < 1e-9)
                        {
                          ixMatch = ix;
                          break;
                        }
                      }

                      if (ixMatch > 0) selTruthXBins.push_back(ixMatch);
                    }

                    const int nyTruth = hTruthTemplate2D->GetNbinsY();
                    for (const int ix : selTruthXBins)
                    {
                      for (int iy = 1; iy <= nyTruth; ++iy)
                      {
                        selTruthGlobalBins.push_back(hTruthTemplate2D->GetBin(ix, iy));
                      }
                    }
                  }

                  if (selTruthXBins.empty() || selTruthGlobalBins.empty()) return s;

                  TH1* hPrev = CloneTH1(hTruthTemplateGlob, "hXJIter0Baseline");
                  if (!hPrev) return s;

                  hPrev->SetDirectory(nullptr);
                  EnsureSumw2(hPrev);
                  hPrev->Reset("ICES");

                  const int nbPrevFill = std::min(hPrev->GetNbinsX(), hMeasuredGlob->GetNbinsX());
                  for (int ib = 1; ib <= nbPrevFill; ++ib)
                  {
                    hPrev->SetBinContent(ib, hMeasuredGlob->GetBinContent(ib));
                    hPrev->SetBinError  (ib, hMeasuredGlob->GetBinError  (ib));
                  }

                  for (int it = 1; it <= kMaxIt; ++it)
                  {
                    RooUnfoldBayes u(&resp, hMeasuredGlob, it);
                    u.SetVerbose(0);
                    u.SetNToys(nToys);

                    TH1* hCurr = nullptr;
                    if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
                    hCurr = u.Hreco(RooUnfold::kCovToy);
                    const auto cov = u.Ereco(RooUnfold::kCovToy);
                    if (gSystem) gSystem->RedirectOutput(0);
                    if (!hCurr) continue;

                    hCurr->SetDirectory(nullptr);
                    EnsureSumw2(hCurr);

                    if (hCurr->GetNbinsX() != hTruthTemplateGlob->GetNbinsX())
                    {
                      delete hCurr;
                      continue;
                    }

                    TH2* hCurr2D = UnflattenGlobalToTH2(
                      hCurr,
                      hTruthTemplate2D,
                      TString::Format("hXJIterCurr2D_%s_it%d", nameSuffix.c_str(), it).Data()
                    );
                    TH2* hPrev2D = UnflattenGlobalToTH2(
                      hPrev,
                      hTruthTemplate2D,
                      TString::Format("hXJIterPrev2D_%s_it%d", nameSuffix.c_str(), it).Data()
                    );

                    if (!hCurr2D || !hPrev2D)
                    {
                      if (hCurr2D) delete hCurr2D;
                      if (hPrev2D) delete hPrev2D;
                      delete hPrev;
                      hPrev = hCurr;
                      continue;
                    }

                    hCurr2D->SetDirectory(nullptr);
                    hPrev2D->SetDirectory(nullptr);

                    double sumV  = 0.0;
                    double sumV2 = 0.0;
                    double sumD2 = 0.0;

                    const int nyTruth = hCurr2D->GetNbinsY();

                    for (const int ix : selTruthXBins)
                    {
                      if (ix < 1 || ix > hCurr2D->GetNbinsX()) continue;

                      for (int iy = 1; iy <= nyTruth; ++iy)
                      {
                        const double v  = hCurr2D->GetBinContent(ix, iy);
                        const double vp = hPrev2D->GetBinContent(ix, iy);
                        const double d  = v - vp;

                        if (v == 0.0 && vp == 0.0) continue;

                        sumV  += v;
                        sumV2 += v * v;
                        sumD2 += d * d;
                      }
                    }

                    double sumCov = 0.0;
                    const int nCovRows = cov.GetNrows();
                    const int nCovCols = cov.GetNcols();

                    for (size_t i = 0; i < selTruthGlobalBins.size(); ++i)
                    {
                      const int ii = selTruthGlobalBins[i];
                      if (ii < 0 || ii >= nCovRows) continue;

                      for (size_t j = 0; j < selTruthGlobalBins.size(); ++j)
                      {
                        const int jj = selTruthGlobalBins[j];
                        if (jj < 0 || jj >= nCovCols) continue;
                        sumCov += cov(ii, jj);
                      }
                    }

                    const double relStat = (sumV > 0.0) ? (std::sqrt(std::max(0.0, sumCov)) / sumV) : 0.0;
                    const double relChg  = (sumV2 > 0.0) ? std::sqrt(sumD2 / sumV2) : 0.0;

                    s.xIt.push_back((double)it);
                    s.exIt.push_back(0.0);
                    s.yRelStat.push_back(relStat);
                    s.eyRelStat.push_back(0.0);
                    s.yRelChange.push_back(relChg);
                    s.eyRelChange.push_back(0.0);
                    s.yQuad.push_back(std::sqrt(relStat * relStat + relChg * relChg));
                    s.eyQuad.push_back(0.0);

                    if (s.yQuad.back() < s.bestQuad)
                    {
                      s.bestQuad = s.yQuad.back();
                      s.bestIt = it;
                    }

                    delete hCurr2D;
                    delete hPrev2D;
                    delete hPrev;
                    hPrev = hCurr;
                  }

                  if (hPrev) delete hPrev;

                  return s;
              };
      
              cout << "\n  [5I] Photon unfolding inputs (photons_INCLUSIVE_QA / PPG12 object-match):\n";
              PrintGet(dsData, phoRecoName,   hPhoRecoData_in);
              PrintGet(dsSim,  phoRecoName,   hPhoRecoSim_in);
              PrintGet(dsSim,  phoTruthName,  hPhoTruthSim_in);
              PrintGet(dsSim,  phoRespName,   hPhoRespSim_in);
              PrintGet(dsSim,  phoFakesName,  hPhoRecoFakesSim_in);
              PrintGet(dsSim,  phoMissesName, hPhoTruthMissesSim_in);
      
              if (!hPhoRecoData_in || !hPhoRecoSim_in || !hPhoTruthSim_in || !hPhoRespSim_in)
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] Missing one or more photons_INCLUSIVE_QA inputs.\n"
                     << "       Need (DATA) h_unfoldRecoPho_pTgamma_ppg12obj and (SIM) h_unfoldRecoPho_pTgamma_ppg12obj, h_unfoldTruthPho_pTgamma_ppg12obj, h2_unfoldResponsePho_pTgamma_ppg12obj.\n"
                     << "       Skipping photons_INCLUSIVE_QA block."
                     << ANSI_RESET << "\n";
                break;
              }
      
              TH1* hPhoRecoData  = CloneTH1(hPhoRecoData_in,  "hPhoRecoData");
              TH1* hPhoRecoSim   = CloneTH1(hPhoRecoSim_in,   "hPhoRecoSim");
              TH1* hPhoTruthSim  = CloneTH1(hPhoTruthSim_in,  "hPhoTruthSim");
              TH2* hPhoRespSim   = CloneTH2(hPhoRespSim_in,   "hPhoRespSim_truthVsReco");
      
              EnsureSumw2(hPhoRecoData);
              EnsureSumw2(hPhoRecoSim);
              EnsureSumw2(hPhoTruthSim);
              EnsureSumw2(hPhoRespSim);
      
                if (gApplyPurityCorrectionForUnfolding)
                {
                  cout << "  [BUILD] DATA :: purity-corrected photon reco input from ABCD counts\n";
                  if (!ApplyPurityCorrectionToRecoPhotonHist(hPhoRecoData))
                  {
                    cout << ANSI_BOLD_YEL
                         << "[WARN] Failed to build purity-corrected photons_INCLUSIVE_QA reco input. Skipping photons_INCLUSIVE_QA block."
                         << ANSI_RESET << "\n";
                    delete hPhoRecoData;
                    delete hPhoRecoSim;
                    delete hPhoTruthSim;
                    delete hPhoRespSim;
                    break;
                  }
                }

                if (gApplyPurityCorrectionForUnfolding)
                {
                  MakePhotonFigure30Like(
                    hPhoRecoData,
                    hPhoRecoSim,
                    phoYieldDir,
                    "ppg12obj",
                    "PPG12 object-match definition"
                  );
                }
        
              TH2D* hPhoResp_measXtruth = TransposeTH2(
                hPhoRespSim,
                "h2_unfoldResponsePho_pTgamma_recoVsTruth",
                "Photon response matrix; p_{T}^{#gamma, reco} [GeV]; p_{T}^{#gamma, truth} [GeV]"
              );
      
            if (!hPhoResp_measXtruth)
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Failed to transpose photons_INCLUSIVE_QA response matrix. Skipping photons_INCLUSIVE_QA block."
                   << ANSI_RESET << "\n";
              delete hPhoRecoData;
              delete hPhoRecoSim;
              delete hPhoTruthSim;
              delete hPhoRespSim;
              break;
            }
      
              const int kDefaultBayesIterPho = 3;
              const int kMaxBayesIterPhoScan = 10;
      
              // Toy settings:
              //   - "final" is for the baseline unfolded spectrum that feeds your main outputs
              //   - "scan" is for iteration-stability / closure / half-closure diagnostics
              const int kNToysPhoFinal = 600;
              const int kNToysPhoScan  = 120;
      
              RooUnfoldResponse respPho(hPhoRecoSim, hPhoTruthSim, hPhoResp_measXtruth, "respPho", "respPho");
      
              const IterScanSummary phoIterScan = BuildPhotonIterScan(
                respPho,
                hPhoRecoData,
                hPhoTruthSim,
                kMaxBayesIterPhoScan,
                kNToysPhoScan
              );
      
              const int kBayesIterPho = (phoIterScan.bestIt > 0 ? phoIterScan.bestIt : kDefaultBayesIterPho);
      
              cout << ANSI_BOLD_CYN
                   << "[PHO ITER AUTO] Using Bayes iteration " << kBayesIterPho;
              if (phoIterScan.bestIt > 0 && std::isfinite(phoIterScan.bestQuad))
              {
                cout << " from quadrature-sum minimum (" << phoIterScan.bestQuad << ")";
              }
              else
              {
                cout << " (fallback default; automatic scan unavailable)";
              }
              cout << ANSI_RESET << "\n";
      
                RooUnfoldBayes    unfoldPhoToy(&respPho, hPhoRecoData, kBayesIterPho);
                unfoldPhoToy.SetVerbose(0);
                unfoldPhoToy.SetNToys(kNToysPhoFinal);
        
                TH1* hPhoUnfoldTruth = nullptr;
                if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
                hPhoUnfoldTruth = unfoldPhoToy.Hreco(RooUnfold::kCovToy);
                if (gSystem) gSystem->RedirectOutput(0);
                if (hPhoUnfoldTruth) hPhoUnfoldTruth->SetDirectory(nullptr);
        
                RooUnfoldBayes    unfoldPhoCov(&respPho, hPhoRecoData, kBayesIterPho);
                unfoldPhoCov.SetVerbose(0);
        
                TH1* hPhoUnfoldTruth_cov = unfoldPhoCov.Hreco(RooUnfold::kCovariance);
                if (hPhoUnfoldTruth_cov) hPhoUnfoldTruth_cov->SetDirectory(nullptr);
      
              // Photon QA outputs
              {
                {
                  TCanvas c("c_pho_resp","c_pho_resp", 900, 750);
                  ApplyCanvasMargins2D(c);
                  c.SetLogz();
      
                  hPhoRespSim->SetTitle("");
                  hPhoRespSim->GetXaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
                  hPhoRespSim->GetYaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                  hPhoRespSim->Draw("colz");
      
                  DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
                  DrawLatexLines(0.14,0.78, { "SIM photon response", "truth #rightarrow reco" }, 0.030, 0.040);
      
                  SaveCanvas(c, JoinPath(phoDir, "pho_response_truthVsReco.png"));
                }
                if (hPhoResp_measXtruth)
                {
                    TCanvas c("c_pho_respT","c_pho_respT", 900, 750);
                    ApplyCanvasMargins2D(c);
                    c.SetLogz();
      
                    hPhoResp_measXtruth->SetTitle("");
                    hPhoResp_measXtruth->GetXaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                    hPhoResp_measXtruth->GetYaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
      
                    const auto& recoEdges  = kUnfoldRecoPtEdges;
                    const auto& truthEdges = kUnfoldTruthPtEdges;
      
                    const double recoMin     = (!recoEdges.empty() ? recoEdges.front() : 8.0);
                    const double recoUFHi    = (recoEdges.size() >= 2 ? recoEdges[1] : 10.0);
                    const double recoOFLo    = (recoEdges.size() >= 2 ? recoEdges[recoEdges.size() - 2] : 35.0);
                    const double recoMax     = (!recoEdges.empty() ? recoEdges.back() : 40.0);
                    const double recoDrawMax = recoOFLo;
      
                    const double truthMin   = (!truthEdges.empty() ? truthEdges.front() : 5.0);
                    const double truthUFMid = (truthEdges.size() >= 2 ? truthEdges[1] : 8.0);
                    const double truthUFHi  = (truthEdges.size() >= 3 ? truthEdges[2] : 10.0);
                    const double truthOFLo  = (truthEdges.size() >= 2 ? truthEdges[truthEdges.size() - 2] : 35.0);
                    const double truthMax   = (!truthEdges.empty() ? truthEdges.back() : 40.0);
      
                    // Show the full truth range, but cut the reco drawing range at 35 GeV
                    // so the reco overflow support bin (35-40) is indicated by the dashed line
                    // at 35 GeV rather than displayed as a full x-axis bin.
                    hPhoResp_measXtruth->GetXaxis()->SetRangeUser(recoMin, recoDrawMax);
                    hPhoResp_measXtruth->GetYaxis()->SetRangeUser(truthMin, truthMax);
      
                    // Debug: verify the response histogram actually contains the expected bin edges
                    cout << "  [pho_response_recoVsTruth] X(reco) nbins=" << hPhoResp_measXtruth->GetXaxis()->GetNbins()
                           << "  firstLowEdge=" << hPhoResp_measXtruth->GetXaxis()->GetBinLowEdge(1)
                           << "  lastUpEdge="   << hPhoResp_measXtruth->GetXaxis()->GetBinUpEdge(hPhoResp_measXtruth->GetXaxis()->GetNbins())
                           << "\n";
                    cout << "  [pho_response_recoVsTruth] Y(truth) nbins=" << hPhoResp_measXtruth->GetYaxis()->GetNbins()
                           << "  firstLowEdge=" << hPhoResp_measXtruth->GetYaxis()->GetBinLowEdge(1)
                           << "  lastUpEdge="   << hPhoResp_measXtruth->GetYaxis()->GetBinUpEdge(hPhoResp_measXtruth->GetYaxis()->GetNbins())
                           << "\n";
      
                    hPhoResp_measXtruth->Draw("colz");
                    if (gPad) gPad->Update();
      
                    TLine lRecoUF(recoUFHi, truthMin, recoUFHi, truthMax);
                    TLine lRecoOF(recoOFLo, truthMin, recoOFLo, truthMax);
                    TLine lTruthUF1(recoMin, truthUFMid, recoDrawMax, truthUFMid);
                    TLine lTruthUF2(recoMin, truthUFHi, recoDrawMax, truthUFHi);
                    TLine lTruthOF(recoMin, truthOFLo, recoDrawMax, truthOFLo);
      
                    lRecoUF.SetLineColor(kBlack);
                    lRecoOF.SetLineColor(kBlack);
                    lTruthUF1.SetLineColor(kBlack);
                    lTruthUF2.SetLineColor(kBlack);
                    lTruthOF.SetLineColor(kBlack);
      
                    lRecoUF.SetLineStyle(2);
                    lRecoOF.SetLineStyle(2);
                    lTruthUF1.SetLineStyle(2);
                    lTruthUF2.SetLineStyle(2);
                    lTruthOF.SetLineStyle(2);
      
                    lRecoUF.SetLineWidth(2);
                    lRecoOF.SetLineWidth(2);
                    lTruthUF1.SetLineWidth(2);
                    lTruthUF2.SetLineWidth(2);
                    lTruthOF.SetLineWidth(2);
      
                    if (recoUFHi > recoMin && recoUFHi < recoMax) lRecoUF.Draw("same");
                    if (recoOFLo > recoMin && recoOFLo < recoMax) lRecoOF.Draw("same");
                    if (truthUFMid > truthMin && truthUFMid < truthMax) lTruthUF1.Draw("same");
                    if (truthUFHi > truthMin && truthUFHi < truthMax) lTruthUF2.Draw("same");
                    if (truthOFLo > truthMin && truthOFLo < truthMax) lTruthOF.Draw("same");
      
                    DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
                    DrawLatexLines(0.14,0.88, { "SIM photon response (transpose)", "reco #rightarrow truth axis order" }, 0.030, 0.040);
                    DrawLatexLines(0.52,0.36,
                                   {
                                     "UF / OF support bins:",
                                     "reco: 8-10 = UF, 35-40 = OF",
                                     "truth: 5-8 and 8-10 = UF, 35-40 = OF"
                                   },
                                   0.026, 0.035);
      
                    c.Modified();
                    c.Update();
      
                    SaveCanvas(c, JoinPath(phoDir, "pho_response_recoVsTruth.png"));
                  }
      
                if (hPhoUnfoldTruth)
                {
                  TH1* hRecoShape = CloneTH1(hPhoRecoData, "hPhoRecoData_forOverlay");
                  TH1* hUnfShape  = CloneTH1(hPhoUnfoldTruth, "hPhoUnfoldTruth_forOverlay");
                  if (hRecoShape && hUnfShape)
                  {
                    hRecoShape->SetDirectory(nullptr);
                    hUnfShape->SetDirectory(nullptr);
      
                    hRecoShape->SetLineColor(2);
                    hRecoShape->SetMarkerColor(2);
                    hRecoShape->SetMarkerStyle(20);
                    hRecoShape->SetLineWidth(2);
      
                    hUnfShape->SetLineColor(1);
                    hUnfShape->SetMarkerColor(1);
                    hUnfShape->SetMarkerStyle(24);
                    hUnfShape->SetLineWidth(2);
      
                    TCanvas c("c_pho_unf","c_pho_unf", 900, 700);
                    ApplyCanvasMargins1D(c);
      
                    const double maxv = std::max(hRecoShape->GetMaximum(), hUnfShape->GetMaximum());
                    hRecoShape->SetMaximum(maxv * 1.35);
                    hRecoShape->SetTitle("");
                    hRecoShape->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    hRecoShape->GetYaxis()->SetTitle("Counts");
      
                    hRecoShape->Draw("E1");
                    hUnfShape->Draw("E1 same");
      
                    TLegend leg(0.55,0.76,0.92,0.90);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.032);
                    leg.AddEntry(hRecoShape, "Run24pp Reco", "lep");
                    leg.AddEntry(hUnfShape,  TString::Format("Unfolded (truth), Bayes it=%d", kBayesIterPho).Data(), "lep");
                    leg.Draw();
      
                    DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsData), 0.034, 0.045);
                    DrawLatexLines(0.14,0.78, { "Photon unfolding: N_{#gamma}(p_{T}^{#gamma})" }, 0.030, 0.040);
      
                    SaveCanvas(c, JoinPath(phoDir, "pho_unfolded_truth_pTgamma_overlay.png"));
      
                    // Also save a log-y version (keep the linear-y output above as-is),
                    // but only draw the common analysis bins (10-35 GeV) and add a ratio subpanel.
                    {
                          const double xPlotMin = 10.0;
                          const double xPlotMax = 35.0;
      
                          auto scanMaxInRange = [&](TH1* h)->double
                          {
                            double maxY = 0.0;
                            if (!h) return maxY;
      
                            const int nb = h->GetNbinsX();
                            for (int ib = 1; ib <= nb; ++ib)
                            {
                              const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                              const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                              if (lo < xPlotMin || hi > xPlotMax) continue;
      
                              const double y  = h->GetBinContent(ib);
                              const double ey = h->GetBinError(ib);
                              if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                              if (y <= 0.0 && ey <= 0.0) continue;
      
                              const double v = y + ey;
                              if (v > maxY) maxY = v;
                            }
                            return maxY;
                          };
      
                          auto scanMinPosInRange = [&](TH1* h)->double
                          {
                            double minY = 1e99;
                            if (!h) return minY;
      
                            const int nb = h->GetNbinsX();
                            for (int ib = 1; ib <= nb; ++ib)
                            {
                              const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                              const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                              if (lo < xPlotMin || hi > xPlotMax) continue;
      
                              const double y = h->GetBinContent(ib);
                              if (!std::isfinite(y)) continue;
                              if (y > 0.0 && y < minY) minY = y;
                            }
                            return minY;
                          };
      
                          auto BuildPhotonAfterOverBeforeRatio = [&](TH1* hAfter, TH1* hBefore)->TH1D*
                          {
                            if (!hAfter || !hBefore) return nullptr;
      
                            vector<double> ratioEdges =
                            {
                              10.0, 12.0, 14.0, 16.0, 18.0,
                              20.0, 22.0, 24.0, 26.0, 35.0
                            };
      
                            TH1D* hRatio = new TH1D(
                              "hPho_afterOverBefore_ratio",
                              "",
                              (int)ratioEdges.size() - 1,
                              &ratioEdges[0]
                            );
                            hRatio->SetDirectory(nullptr);
                            hRatio->Sumw2();
      
                            for (int ib = 1; ib <= hRatio->GetNbinsX(); ++ib)
                            {
                              const double lo  = hRatio->GetXaxis()->GetBinLowEdge(ib);
                              const double hi  = hRatio->GetXaxis()->GetBinUpEdge(ib);
                              const double cen = hRatio->GetXaxis()->GetBinCenter(ib);
      
                              const int iAfter  = hAfter ->GetXaxis()->FindBin(cen);
                              const int iBefore = hBefore->GetXaxis()->FindBin(cen);
      
                              if (iAfter  < 1 || iAfter  > hAfter ->GetNbinsX()) continue;
                              if (iBefore < 1 || iBefore > hBefore->GetNbinsX()) continue;
      
                              const double afterLo  = hAfter ->GetXaxis()->GetBinLowEdge(iAfter);
                              const double afterHi  = hAfter ->GetXaxis()->GetBinUpEdge(iAfter);
                              const double beforeLo = hBefore->GetXaxis()->GetBinLowEdge(iBefore);
                              const double beforeHi = hBefore->GetXaxis()->GetBinUpEdge(iBefore);
      
                              if (std::fabs(afterLo  - lo) > 1e-6 || std::fabs(afterHi  - hi) > 1e-6) continue;
                              if (std::fabs(beforeLo - lo) > 1e-6 || std::fabs(beforeHi - hi) > 1e-6) continue;
      
                              const double num  = hAfter ->GetBinContent(iAfter);
                              const double eNum = hAfter ->GetBinError  (iAfter);
                              const double den  = hBefore->GetBinContent(iBefore);
                              const double eDen = hBefore->GetBinError  (iBefore);
      
                              if (!std::isfinite(num) || !std::isfinite(eNum) ||
                                  !std::isfinite(den) || !std::isfinite(eDen) || den <= 0.0)
                              {
                                continue;
                              }
      
                              const double val = num / den;
                              const double var = (eNum * eNum) / (den * den)
                                               + (num * num * eDen * eDen) / (den * den * den * den);
                              const double err = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : 0.0;
      
                              hRatio->SetBinContent(ib, val);
                              hRatio->SetBinError  (ib, err);
                            }
      
                            return hRatio;
                          };
      
                          const double maxvTop = std::max(scanMaxInRange(hRecoShape), scanMaxInRange(hUnfShape));
      
                          double minPos = std::min(scanMinPosInRange(hRecoShape), scanMinPosInRange(hUnfShape));
                          if (!(minPos < 1e98)) minPos = 1e-3;
      
                          TH1D* hRatio = BuildPhotonAfterOverBeforeRatio(hUnfShape, hRecoShape);
      
                          double ratioMin = 0.8;
                          double ratioMax = 1.2;
                          if (hRatio)
                          {
                            double rMin =  1e99;
                            double rMax = -1e99;
      
                            for (int ib = 1; ib <= hRatio->GetNbinsX(); ++ib)
                            {
                              const double y  = hRatio->GetBinContent(ib);
                              const double ey = hRatio->GetBinError(ib);
                              if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                              if (y <= 0.0 && ey <= 0.0) continue;
      
                              rMin = std::min(rMin, y - ey);
                              rMax = std::max(rMax, y + ey);
                            }
      
                            if (rMin < 1e98 && rMax > -1e98 && rMin < rMax)
                            {
                              const double span = std::max(rMax - rMin, 0.04);
                              const double pad  = 0.18 * span;
                              ratioMin = rMin - pad;
                              ratioMax = rMax + pad;
      
                              if (ratioMin > 1.0) ratioMin = 1.0 - 0.5 * span;
                              if (ratioMax < 1.0) ratioMax = 1.0 + 0.5 * span;
                              if (ratioMin < 0.0) ratioMin = 0.0;
                            }
                          }
      
                          c.Clear();
      
                          TPad* pTop = new TPad("pTop_pho_unf_logy", "pTop_pho_unf_logy", 0.0, 0.30, 1.0, 1.0);
                          TPad* pBot = new TPad("pBot_pho_unf_logy", "pBot_pho_unf_logy", 0.0, 0.00, 1.0, 0.30);
      
                          pTop->SetLeftMargin(0.12);
                          pTop->SetRightMargin(0.04);
                          pTop->SetTopMargin(0.08);
                          pTop->SetBottomMargin(0.02);
                          pTop->SetTicks(1, 1);
                          pTop->SetLogy(1);
      
                          pBot->SetLeftMargin(0.12);
                          pBot->SetRightMargin(0.04);
                          pBot->SetTopMargin(0.02);
                          pBot->SetBottomMargin(0.30);
                          pBot->SetTicks(1, 1);
                          pBot->SetGridy(1);
      
                          pTop->Draw();
                          pBot->Draw();
      
                          pTop->cd();
      
                          hRecoShape->SetTitle("");
                          hRecoShape->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                          hUnfShape->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
      
                          hRecoShape->GetXaxis()->SetTitle("");
                          hRecoShape->GetXaxis()->SetLabelSize(0.0);
                          hRecoShape->GetXaxis()->SetTitleSize(0.0);
      
                          hRecoShape->GetYaxis()->SetTitle("Counts");
                          hRecoShape->GetYaxis()->SetTitleSize(0.055);
                          hRecoShape->GetYaxis()->SetTitleOffset(0.95);
                          hRecoShape->GetYaxis()->SetLabelSize(0.045);
      
                          hRecoShape->SetMinimum(std::max(minPos * 0.5, 1e-6));
                          hRecoShape->SetMaximum((maxvTop > 0.0) ? (maxvTop * 3.0) : 1.0);
      
                          hRecoShape->Draw("E1");
                          hUnfShape->Draw("E1 same");
      
                          TLegend leg(0.55,0.76,0.92,0.90);
                          leg.SetTextFont(42);
                          leg.SetTextSize(0.032);
                          leg.AddEntry(hRecoShape, "PP DATA (reco)", "lep");
                          leg.AddEntry(hUnfShape,  TString::Format("Unfolded (truth), Bayes it=%d", kBayesIterPho).Data(), "lep");
                          leg.Draw();
      
                          DrawLatexLines(0.14,0.92,
                                         { "Photon 1D Unfolding, N_{#gamma}(p_{T}^{#gamma}), Photon 4 GeV + MBD NS #geq 1" },
                                         0.034, 0.045);
      
                          pBot->cd();
      
                          if (hRatio)
                          {
                            hRatio->SetTitle("");
                            hRatio->SetMarkerStyle(20);
                            hRatio->SetMarkerSize(0.95);
                            hRatio->SetLineWidth(2);
                            hRatio->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                            hRatio->GetYaxis()->SetTitle("After / before unfolding");
                            hRatio->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                            hRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
                            hRatio->GetXaxis()->SetTitleSize(0.12);
                            hRatio->GetXaxis()->SetLabelSize(0.11);
                            hRatio->GetXaxis()->SetTitleOffset(1.00);
                            hRatio->GetYaxis()->SetTitleSize(0.10);
                            hRatio->GetYaxis()->SetLabelSize(0.09);
                            hRatio->GetYaxis()->SetTitleOffset(0.55);
                            hRatio->GetYaxis()->SetNdivisions(505);
                            hRatio->Draw("E1");
      
                            TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                            l1.SetLineStyle(2);
                            l1.SetLineWidth(2);
                            l1.Draw("same");
                          }
                          else
                          {
                            TH1F frame("frame_pho_ratio","", 1, xPlotMin, xPlotMax);
                            frame.SetMinimum(ratioMin);
                            frame.SetMaximum(ratioMax);
                            frame.SetTitle("");
                            frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                            frame.GetYaxis()->SetTitle("After / before unfolding");
                            frame.GetXaxis()->SetTitleSize(0.12);
                            frame.GetXaxis()->SetLabelSize(0.11);
                            frame.GetXaxis()->SetTitleOffset(1.00);
                            frame.GetYaxis()->SetTitleSize(0.10);
                            frame.GetYaxis()->SetLabelSize(0.09);
                            frame.GetYaxis()->SetTitleOffset(0.55);
                            frame.GetYaxis()->SetNdivisions(505);
                            frame.Draw("axis");
      
                            TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                            l1.SetLineStyle(2);
                            l1.SetLineWidth(2);
                            l1.Draw("same");
                          }
      
                          c.cd();
                          c.Modified();
                          c.Update();
                          SaveCanvas(c, JoinPath(phoDir, "pho_unfolded_truth_pTgamma_overlay_logy.png"));
      
                          if (hRatio) delete hRatio;
                      }
      
                      // photon efficiency/purity diagnostics
                      //
                      // These are SIM truth-matching bookkeeping diagnostics that explain the
                      // normalization gap in before/after unfolding overlays.
                      //
                      // (1) Purity (SIM, reco space):
                      //     purity(pT) = 1 - N_fakeReco(pT)/N_reco(pT)
                      //
                      // (2) Efficiency (SIM, truth space):
                      //     eff(pT)    = 1 - N_missTruth(pT)/N_truth(pT)
                      //
                      // (3) Truth/Reconstruction scale factor (SIM):
                      //     N_truth(pT)/N_reco(pT)  (with bin-mapping when axes differ)
                      //
                      // (4) Implied "efficiency-like" curve in DATA (not purely data-driven):
                      //     eps_eff,data(pT) = N_reco,data(pT) / N_truth,data(unfolded)(pT)
                      //     (computed with bin-mapping when axes differ)
                      //
                      // Outputs (to <phoDir>, i.e. unfolding/radii/<rXX>/photons):
                      //   - pho_efficiencyEff_data_vs_pTgamma.png
                      //   - pho_purity_sim_vs_pTgamma.png                 (if fakes hist exists)
                      //   - pho_efficiency_sim_vs_pTgamma.png             (if misses hist exists)
                      //   - pho_truthOverReco_sim_vs_pTgamma.png          (bin-mapped)
                      //   - pho_recoOverTruth_sim_vs_pTgamma.png          (bin-mapped)
                      //   - pho_efficiencyEff_data_vs_efficiency_sim.png  (overlay, if both exist)
                      // -------------------------------------------------------------------
                      {
                        auto sameBinning = [&](TH1* a, TH1* b)->bool
                        {
                          if (!a || !b) return false;
                          if (a->GetNbinsX() != b->GetNbinsX()) return false;
                          const int nb = a->GetNbinsX();
                          for (int ib = 1; ib <= nb + 1; ++ib)
                          {
                            const double ea = a->GetXaxis()->GetBinUpEdge(ib);
                            const double eb = b->GetXaxis()->GetBinUpEdge(ib);
                            if (std::fabs(ea - eb) > 1e-9) return false;
                          }
                          return true;
                        };
      
                        auto mapToRefBinning = [&](TH1* src, TH1* ref, const char* newName)->TH1*
                        {
                          if (!src || !ref) return nullptr;
      
                          TH1* h = CloneTH1(ref, newName);
                          if (!h) return nullptr;
      
                          h->SetDirectory(nullptr);
                          EnsureSumw2(h);
                          h->Reset("ICES");
      
                          const int nb = ref->GetNbinsX();
                          int nBad = 0;
      
                          for (int ib = 1; ib <= nb; ++ib)
                          {
                            const double x  = ref->GetXaxis()->GetBinCenter(ib);
                            const int isrc  = src->GetXaxis()->FindBin(x);
      
                            const double xsLo = src->GetXaxis()->GetBinLowEdge(isrc);
                            const double xsHi = src->GetXaxis()->GetBinUpEdge(isrc);
                            const double xrLo = ref->GetXaxis()->GetBinLowEdge(ib);
                            const double xrHi = ref->GetXaxis()->GetBinUpEdge(ib);
      
                            if (std::fabs(xsLo - xrLo) > 1e-3 || std::fabs(xsHi - xrHi) > 1e-3) ++nBad;
      
                            h->SetBinContent(ib, src->GetBinContent(isrc));
                            h->SetBinError  (ib, src->GetBinError  (isrc));
                          }
      
                          if (nBad > 0)
                          {
                            cout << ANSI_BOLD_YEL
                                 << "[WARN] Photon diagnostics: mapped histogram '" << newName
                                 << "' had " << nBad << " bins with mismatched edges (copied by bin-center)."
                                 << ANSI_RESET << "\n";
                          }
      
                          return h;
                        };
      
                        auto drawLineAtOne = [&](TH1* h)
                        {
                          if (!h) return;
                          const double xmin = h->GetXaxis()->GetXmin();
                          const double xmax = h->GetXaxis()->GetXmax();
                          TLine l1(xmin, 1.0, xmax, 1.0);
                          l1.SetLineStyle(2);
                          l1.SetLineWidth(2);
                          l1.Draw("same");
                        };
      
                        // -----------------------------
                        // DATA: eps_eff = N_reco,data / N_truth,data(unfolded)
                        // -----------------------------
                        TH1* hEffData = nullptr;
                        {
                          TH1* hRecoD  = hPhoRecoData;
                          TH1* hTruthD = hPhoUnfoldTruth;
      
                          if (hRecoD && hTruthD)
                          {
                            TH1* hRecoD_m = nullptr;
      
                            if (sameBinning(hRecoD, hTruthD))
                            {
                              hRecoD_m = CloneTH1(hRecoD, "h_pho_recoData_counts_forEff");
                              if (hRecoD_m) { hRecoD_m->SetDirectory(nullptr); EnsureSumw2(hRecoD_m); }
                            }
                            else
                            {
                              hRecoD_m = mapToRefBinning(hRecoD, hTruthD, "h_pho_recoData_counts_mappedToTruthBins_forEff");
                            }
      
                            if (hRecoD_m)
                            {
                              hEffData = CloneTH1(hRecoD_m, "h_pho_effData_recoOverUnfoldTruth");
                              if (hEffData)
                              {
                                hEffData->SetDirectory(nullptr);
                                EnsureSumw2(hEffData);
                                hEffData->Divide(hTruthD);
      
                                hEffData->SetTitle("");
                                hEffData->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                                hEffData->GetYaxis()->SetTitle("#epsilon_{#gamma}^{eff,data} = N_{#gamma}^{reco,data} / N_{#gamma}^{truth,data (unfolded)}");
                                hEffData->SetMarkerStyle(20);
                                hEffData->SetMarkerSize(1.1);
                                hEffData->SetLineWidth(2);
      
                                TCanvas c("c_pho_effData", "c_pho_effData", 900, 700);
                                ApplyCanvasMargins1D(c);
      
                                hEffData->GetYaxis()->SetRangeUser(0.0, 1.2);
                                hEffData->Draw("E1");
                                drawLineAtOne(hEffData);
      
                                DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsData), 0.034, 0.045);
                                DrawLatexLines(0.14,0.78, { "Photon unfolding diagnostic: implied #epsilon_{#gamma}^{eff,data}(p_{T}^{#gamma})" }, 0.030, 0.040);
      
                                SaveCanvas(c, JoinPath(phoDir, "pho_efficiencyEff_data_vs_pTgamma.png"));
                              }
                            }
      
                            if (hRecoD_m && hRecoD_m != hRecoD) delete hRecoD_m;
                          }
                        }
      
                        // -----------------------------
                        // SIM: purity, efficiency, truth/reco ratios
                        // -----------------------------
                        TH1* hPurSim = nullptr;
                        TH1* hEffSim = nullptr;
                        TH1* hTruthOverRecoSim = nullptr;
                        TH1* hRecoOverTruthSim = nullptr;
      
                          // Purity: 1 - fakes/reco  (reco space)
                          if (hPhoRecoSim && hPhoRecoFakesSim_in)
                          {
                            TH1* hFakeOverReco = CloneTH1(hPhoRecoFakesSim_in, "h_pho_fakeOverReco_sim");
                            if (hFakeOverReco)
                            {
                              hFakeOverReco->SetDirectory(nullptr);
                              EnsureSumw2(hFakeOverReco);
                              hFakeOverReco->Divide(hPhoRecoSim);
      
                            hPurSim = CloneTH1(hFakeOverReco, "h_pho_purity_sim");
                            if (hPurSim)
                            {
                              hPurSim->SetDirectory(nullptr);
                              EnsureSumw2(hPurSim);
      
                              for (int ib = 1; ib <= hPurSim->GetNbinsX(); ++ib)
                              {
                                const double v  = 1.0 - hFakeOverReco->GetBinContent(ib);
                                const double ev = hFakeOverReco->GetBinError(ib);
                                hPurSim->SetBinContent(ib, v);
                                hPurSim->SetBinError  (ib, ev);
                              }
      
                              hPurSim->SetTitle("");
                              hPurSim->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                              hPurSim->GetYaxis()->SetTitle("Purity(p_{T}^{#gamma}) = 1 - N_{fake}^{reco}/N_{#gamma}^{reco}");
                              hPurSim->SetMarkerStyle(21);
                              hPurSim->SetMarkerSize(1.1);
                              hPurSim->SetLineWidth(2);
      
                              TCanvas c("c_pho_purity_sim", "c_pho_purity_sim", 900, 700);
                              ApplyCanvasMargins1D(c);
      
                              hPurSim->GetYaxis()->SetRangeUser(0.0, 1.2);
                              hPurSim->Draw("E1");
                              drawLineAtOne(hPurSim);
      
                              DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
                              DrawLatexLines(0.14,0.78, { "SIM photon purity: 1 - (reco fakes)/(reco selected)" }, 0.030, 0.040);
      
                              SaveCanvas(c, JoinPath(phoDir, "pho_purity_sim_vs_pTgamma.png"));
                            }
      
                            delete hFakeOverReco;
                          }
                        }
                        else
                        {
                          cout << ANSI_BOLD_YEL
                               << "[WARN] Photon purity plot: missing SIM fakes histogram (h_unfoldRecoPhoFakes_pTgamma). Skipping purity plot."
                               << ANSI_RESET << "\n";
                        }
      
                          // Efficiency: 1 - misses/truth (truth space)
                        if (hPhoTruthSim && hPhoTruthMissesSim_in)
                        {
                            TH1* hMissOverTruth = CloneTH1(hPhoTruthMissesSim_in, "h_pho_missOverTruth_sim");
                            if (hMissOverTruth)
                            {
                              hMissOverTruth->SetDirectory(nullptr);
                              EnsureSumw2(hMissOverTruth);
                              hMissOverTruth->Divide(hPhoTruthSim);
      
                            hEffSim = CloneTH1(hMissOverTruth, "h_pho_efficiency_sim");
                            if (hEffSim)
                            {
                              hEffSim->SetDirectory(nullptr);
                              EnsureSumw2(hEffSim);
      
                              for (int ib = 1; ib <= hEffSim->GetNbinsX(); ++ib)
                              {
                                const double v  = 1.0 - hMissOverTruth->GetBinContent(ib);
                                const double ev = hMissOverTruth->GetBinError(ib);
                                hEffSim->SetBinContent(ib, v);
                                hEffSim->SetBinError  (ib, ev);
                              }
      
                              hEffSim->SetTitle("");
                              hEffSim->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                              hEffSim->GetYaxis()->SetTitle("#epsilon_{#gamma}^{MC}(p_{T}^{#gamma}) = 1 - N_{miss}^{truth}/N_{#gamma}^{truth}");
                              hEffSim->SetMarkerStyle(22);
                              hEffSim->SetMarkerSize(1.1);
                              hEffSim->SetLineWidth(2);
      
                              TCanvas c("c_pho_efficiency_sim", "c_pho_efficiency_sim", 900, 700);
                              ApplyCanvasMargins1D(c);
      
                              hEffSim->GetYaxis()->SetRangeUser(0.0, 1.2);
                              hEffSim->Draw("E1");
                              drawLineAtOne(hEffSim);
      
                              DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
                              DrawLatexLines(0.14,0.78, { "SIM photon efficiency: 1 - (truth misses)/(truth signal)" }, 0.030, 0.040);
      
                              SaveCanvas(c, JoinPath(phoDir, "pho_efficiency_sim_vs_pTgamma.png"));
                            }
      
                            delete hMissOverTruth;
                          }
                        }
                        else
                        {
                          cout << ANSI_BOLD_YEL
                               << "[WARN] Photon efficiency plot: missing SIM misses histogram (h_unfoldTruthPhoMisses_pTgamma). Skipping efficiency plot."
                               << ANSI_RESET << "\n";
                        }
      
                        // Truth/Reco ratios in SIM (bin-mapped so it is well-defined)
                        if (hPhoTruthSim && hPhoRecoSim)
                        {
                          TH1* hRecoSim_mTruth = nullptr;
                          if (sameBinning(hPhoRecoSim, hPhoTruthSim))
                          {
                            hRecoSim_mTruth = CloneTH1(hPhoRecoSim, "h_pho_recoSim_counts_forTruthOverReco");
                            if (hRecoSim_mTruth) { hRecoSim_mTruth->SetDirectory(nullptr); EnsureSumw2(hRecoSim_mTruth); }
                          }
                          else
                          {
                            hRecoSim_mTruth = mapToRefBinning(hPhoRecoSim, hPhoTruthSim, "h_pho_recoSim_counts_mappedToTruthBins_forTruthOverReco");
                          }
      
                          if (hRecoSim_mTruth)
                          {
                            hTruthOverRecoSim = CloneTH1(hPhoTruthSim, "h_pho_truthOverReco_sim");
                            if (hTruthOverRecoSim)
                            {
                              hTruthOverRecoSim->SetDirectory(nullptr);
                              EnsureSumw2(hTruthOverRecoSim);
                              hTruthOverRecoSim->Divide(hRecoSim_mTruth);
      
                              hTruthOverRecoSim->SetTitle("");
                              hTruthOverRecoSim->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                              hTruthOverRecoSim->GetYaxis()->SetTitle("N_{#gamma}^{truth} / N_{#gamma}^{reco}  (SIM, bin-mapped)");
                              hTruthOverRecoSim->SetMarkerStyle(20);
                              hTruthOverRecoSim->SetMarkerSize(1.1);
                              hTruthOverRecoSim->SetLineWidth(2);
      
                              TCanvas c("c_pho_truthOverReco_sim", "c_pho_truthOverReco_sim", 900, 700);
                              ApplyCanvasMargins1D(c);
      
                              hTruthOverRecoSim->Draw("E1");
      
                              DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
                              DrawLatexLines(0.14,0.78, { "SIM normalization lever-arm: N_{#gamma}^{truth} / N_{#gamma}^{reco}" }, 0.030, 0.040);
      
                              SaveCanvas(c, JoinPath(phoDir, "pho_truthOverReco_sim_vs_pTgamma.png"));
                            }
      
                            if (hRecoSim_mTruth && hRecoSim_mTruth != hPhoRecoSim) delete hRecoSim_mTruth;
                          }
      
                          TH1* hTruthSim_mReco = nullptr;
                          if (sameBinning(hPhoTruthSim, hPhoRecoSim))
                          {
                            hTruthSim_mReco = CloneTH1(hPhoTruthSim, "h_pho_truthSim_counts_forRecoOverTruth");
                            if (hTruthSim_mReco) { hTruthSim_mReco->SetDirectory(nullptr); EnsureSumw2(hTruthSim_mReco); }
                          }
                          else
                          {
                            hTruthSim_mReco = mapToRefBinning(hPhoTruthSim, hPhoRecoSim, "h_pho_truthSim_counts_mappedToRecoBins_forRecoOverTruth");
                          }
      
                          if (hTruthSim_mReco)
                          {
                            hRecoOverTruthSim = CloneTH1(hPhoRecoSim, "h_pho_recoOverTruth_sim");
                            if (hRecoOverTruthSim)
                            {
                              hRecoOverTruthSim->SetDirectory(nullptr);
                              EnsureSumw2(hRecoOverTruthSim);
                              hRecoOverTruthSim->Divide(hTruthSim_mReco);
      
                              hRecoOverTruthSim->SetTitle("");
                              hRecoOverTruthSim->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                              hRecoOverTruthSim->GetYaxis()->SetTitle("N_{#gamma}^{reco} / N_{#gamma}^{truth}  (SIM, bin-mapped)");
                              hRecoOverTruthSim->SetMarkerStyle(20);
                              hRecoOverTruthSim->SetMarkerSize(1.1);
                              hRecoOverTruthSim->SetLineWidth(2);
      
                              TCanvas c("c_pho_recoOverTruth_sim", "c_pho_recoOverTruth_sim", 900, 700);
                              ApplyCanvasMargins1D(c);
      
                              hRecoOverTruthSim->GetYaxis()->SetRangeUser(0.0, 1.2);
                              hRecoOverTruthSim->Draw("E1");
                              drawLineAtOne(hRecoOverTruthSim);
      
                              DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
                              DrawLatexLines(0.14,0.78, { "SIM efficiency-like: N_{#gamma}^{reco} / N_{#gamma}^{truth}" }, 0.030, 0.040);
      
                              SaveCanvas(c, JoinPath(phoDir, "pho_recoOverTruth_sim_vs_pTgamma.png"));
                            }
      
                            if (hTruthSim_mReco && hTruthSim_mReco != hPhoTruthSim) delete hTruthSim_mReco;
                          }
                        }
      
                        // -----------------------------
                        // Overlay: eps_eff,data vs eps_MC (if both exist)
                        // -----------------------------
                        if (hEffData && hEffSim)
                        {
                          TCanvas c("c_pho_effData_vs_effSim", "c_pho_effData_vs_effSim", 900, 700);
                          ApplyCanvasMargins1D(c);
      
                          TH1* hFrame = CloneTH1(hEffSim, "hFrame_effData_vs_effSim");
                          if (hFrame)
                          {
                            hFrame->SetDirectory(nullptr);
                            hFrame->Reset("ICES");
                            hFrame->SetTitle("");
                            hFrame->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                            hFrame->GetYaxis()->SetTitle("#epsilon_{#gamma}(p_{T}^{#gamma})");
                            hFrame->GetYaxis()->SetRangeUser(0.0, 1.2);
                            hFrame->Draw("axis");
                          }
      
                          hEffSim->SetMarkerColor(kBlue + 1);
                          hEffSim->SetLineColor(kBlue + 1);
                          hEffData->SetMarkerColor(kBlack);
                          hEffData->SetLineColor(kBlack);
      
                          hEffSim->Draw("E1 same");
                          hEffData->Draw("E1 same");
      
                          TGraphErrors gLegData(1), gLegSim(1);
                          {
                            int ib = 1;
                            const double x  = hEffData->GetXaxis()->GetBinCenter(ib);
                            const double y  = hEffData->GetBinContent(ib);
                            const double ey = hEffData->GetBinError(ib);
                            gLegData.SetPoint(0, x, y);
                            gLegData.SetPointError(0, 0.0, ey);
                            gLegData.SetMarkerStyle(hEffData->GetMarkerStyle());
                            gLegData.SetMarkerSize(hEffData->GetMarkerSize());
                            gLegData.SetMarkerColor(hEffData->GetMarkerColor());
                            gLegData.SetLineColor(hEffData->GetLineColor());
                            gLegData.SetLineWidth(hEffData->GetLineWidth());
                          }
                          {
                            int ib = 1;
                            const double x  = hEffSim->GetXaxis()->GetBinCenter(ib);
                            const double y  = hEffSim->GetBinContent(ib);
                            const double ey = hEffSim->GetBinError(ib);
                            gLegSim.SetPoint(0, x, y);
                            gLegSim.SetPointError(0, 0.0, ey);
                            gLegSim.SetMarkerStyle(hEffSim->GetMarkerStyle());
                            gLegSim.SetMarkerSize(hEffSim->GetMarkerSize());
                            gLegSim.SetMarkerColor(hEffSim->GetMarkerColor());
                            gLegSim.SetLineColor(hEffSim->GetLineColor());
                            gLegSim.SetLineWidth(hEffSim->GetLineWidth());
                          }
      
                            TLegend leg(0.24, 0.74, 0.5, 0.88);
                            leg.SetBorderSize(0);
                            leg.SetFillStyle(0);
                            leg.SetTextFont(42);
                            leg.SetTextSize(0.04);
                            leg.AddEntry(&gLegData, "#epsilon_{#gamma}^{eff,data} (reco / unfolded truth)", "pe");
                            leg.AddEntry(&gLegSim,  "#epsilon_{#gamma}^{MC} (1 - misses/truth)", "pe");
                            leg.Draw();
        
                            {
                              const double xmin = hFrame ? hFrame->GetXaxis()->GetXmin() : hEffSim->GetXaxis()->GetXmin();
                              const double xmax = hFrame ? hFrame->GetXaxis()->GetXmax() : hEffSim->GetXaxis()->GetXmax();
                              TLine l1(xmin, 1.0, xmax, 1.0);
                              l1.SetLineStyle(2);
                              l1.SetLineWidth(2);
                              l1.SetLineColor(kRed + 1);
                              l1.Draw("same");
                            }
        
                            {
                              TLatex tx;
                              tx.SetNDC();
                              tx.SetTextFont(42);
                              tx.SetTextAlign(22);
                              tx.SetTextSize(0.040);
                              tx.DrawLatex(0.50, 0.965, "Photon efficiency diagnostic: DATA reco / unfolded truth vs SIM (1 - misses / truth)");
                            }
        
                            SaveCanvas(c, JoinPath(phoDir, "pho_efficiencyEff_data_vs_efficiency_sim.png"));
      
                          if (hFrame) delete hFrame;
                        }
      
                        if (hRecoOverTruthSim) delete hRecoOverTruthSim;
                        if (hTruthOverRecoSim) delete hTruthOverRecoSim;
                        if (hEffSim) delete hEffSim;
                        if (hPurSim) delete hPurSim;
                        if (hEffData) delete hEffData;
                      }
      
                      delete hRecoShape;
                      delete hUnfShape;
                  }
                }
              }
              // Photon QA outputs
              {
                auto sameBinning = [&](TH1* a, TH1* b)->bool
                {
                  if (!a || !b) return false;
                  if (a->GetNbinsX() != b->GetNbinsX()) return false;
                  const int nb = a->GetNbinsX();
                  for (int ib = 1; ib <= nb + 1; ++ib)
                  {
                    const double ea = a->GetXaxis()->GetBinUpEdge(ib);
                    const double eb = b->GetXaxis()->GetBinUpEdge(ib);
                    if (std::fabs(ea - eb) > 1e-9) return false;
                  }
                  return true;
                };
      
                auto mapToRefBinning = [&](TH1* src, TH1* ref, const char* newName)->TH1*
                {
                  if (!src || !ref) return nullptr;
      
                  TH1* h = CloneTH1(ref, newName);
                  if (!h) return nullptr;
      
                  h->SetDirectory(nullptr);
                  EnsureSumw2(h);
                  h->Reset("ICES");
      
                  const int nb = ref->GetNbinsX();
                  int nBad = 0;
      
                  for (int ib = 1; ib <= nb; ++ib)
                  {
                    const double x  = ref->GetXaxis()->GetBinCenter(ib);
                    const int isrc  = src->GetXaxis()->FindBin(x);
      
                    const double xsLo = src->GetXaxis()->GetBinLowEdge(isrc);
                    const double xsHi = src->GetXaxis()->GetBinUpEdge(isrc);
                    const double xrLo = ref->GetXaxis()->GetBinLowEdge(ib);
                    const double xrHi = ref->GetXaxis()->GetBinUpEdge(ib);
      
                    if (std::fabs(xsLo - xrLo) > 1e-3 || std::fabs(xsHi - xrHi) > 1e-3) ++nBad;
      
                    h->SetBinContent(ib, src->GetBinContent(isrc));
                    h->SetBinError  (ib, src->GetBinError  (isrc));
                  }
      
                  if (nBad > 0)
                  {
                    cout << ANSI_BOLD_YEL
                         << "[WARN] Photon QA: mapped histogram '" << newName
                         << "' had " << nBad << " bins with mismatched edges (copied by bin-center)."
                         << ANSI_RESET << "\n";
                  }
      
                  return h;
                };
      
                auto cloneOrMapToRef = [&](TH1* src, TH1* ref, const char* newName)->TH1*
                {
                  if (!src || !ref) return nullptr;
      
                  TH1* h = nullptr;
                  if (sameBinning(src, ref))
                  {
                    h = CloneTH1(src, newName);
                    if (h)
                    {
                      h->SetDirectory(nullptr);
                      EnsureSumw2(h);
                    }
                  }
                  else
                  {
                    h = mapToRefBinning(src, ref, newName);
                  }
                  return h;
                };
      
                auto makeRatioOnRef = [&](TH1* hNum, TH1* hDen, TH1* hRef, const char* newName)->TH1*
                {
                  if (!hNum || !hDen || !hRef) return nullptr;
      
                  TH1* hNumR = cloneOrMapToRef(hNum, hRef, TString::Format("%s_numRef", newName).Data());
                  TH1* hDenR = cloneOrMapToRef(hDen, hRef, TString::Format("%s_denRef", newName).Data());
                  if (!hNumR || !hDenR)
                  {
                    if (hNumR) delete hNumR;
                    if (hDenR) delete hDenR;
                    return nullptr;
                  }
      
                  TH1* hR = CloneTH1(hNumR, newName);
                  if (hR)
                  {
                    hR->SetDirectory(nullptr);
                    EnsureSumw2(hR);
                    hR->Divide(hDenR);
                  }
      
                  delete hNumR;
                  delete hDenR;
                  return hR;
                };
      
                auto makeOneMinus = [&](TH1* hIn, const char* newName)->TH1*
                {
                  if (!hIn) return nullptr;
      
                  TH1* h = CloneTH1(hIn, newName);
                  if (!h) return nullptr;
      
                  h->SetDirectory(nullptr);
                  EnsureSumw2(h);
      
                  for (int ib = 0; ib <= h->GetNbinsX() + 1; ++ib)
                  {
                    if (ib == 0 || ib == h->GetNbinsX() + 1)
                    {
                      h->SetBinContent(ib, 0.0);
                      h->SetBinError  (ib, 0.0);
                      continue;
                    }
      
                    const double v  = 1.0 - hIn->GetBinContent(ib);
                    const double ev = hIn->GetBinError(ib);
                    h->SetBinContent(ib, v);
                    h->SetBinError  (ib, ev);
                  }
      
                  return h;
                };
      
                auto makeHistFromAxis = [&](const TAxis* ax, const char* newName, const char* yTitle)->TH1D*
                {
                  if (!ax) return nullptr;
      
                  TH1D* h = nullptr;
                  if (ax->GetXbins() && ax->GetXbins()->GetSize() > 0)
                  {
                    h = new TH1D(newName, "", ax->GetNbins(), ax->GetXbins()->GetArray());
                  }
                  else
                  {
                    h = new TH1D(newName, "", ax->GetNbins(), ax->GetXmin(), ax->GetXmax());
                  }
      
                  if (!h) return nullptr;
                  h->SetDirectory(nullptr);
                  h->Sumw2();
                  h->SetTitle("");
                  h->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                  h->GetYaxis()->SetTitle(yTitle);
                  return h;
                };
      
                auto normalizeResponseByTruth = [&](TH2* hIn, const char* newName)->TH2*
                {
                  if (!hIn) return nullptr;
      
                  TH2* h = CloneTH2(hIn, newName);
                  if (!h) return nullptr;
      
                  h->SetDirectory(nullptr);
                  EnsureSumw2(h);
                  h->Reset("ICES");
      
                  for (int ix = 1; ix <= hIn->GetNbinsX(); ++ix)
                  {
                    double sum = 0.0;
                    for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy) sum += hIn->GetBinContent(ix, iy);
                    if (!(sum > 0.0)) continue;
      
                    for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy)
                    {
                      const double v = hIn->GetBinContent(ix, iy) / sum;
                      h->SetBinContent(ix, iy, v);
                      h->SetBinError  (ix, iy, 0.0);
                    }
                  }
      
                  return h;
                };
      
                auto normalizeResponseByReco = [&](TH2* hIn, const char* newName)->TH2*
                {
                  if (!hIn) return nullptr;
      
                  TH2* h = CloneTH2(hIn, newName);
                  if (!h) return nullptr;
      
                  h->SetDirectory(nullptr);
                  EnsureSumw2(h);
                  h->Reset("ICES");
      
                  for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy)
                  {
                    double sum = 0.0;
                    for (int ix = 1; ix <= hIn->GetNbinsX(); ++ix) sum += hIn->GetBinContent(ix, iy);
                    if (!(sum > 0.0)) continue;
      
                    for (int ix = 1; ix <= hIn->GetNbinsX(); ++ix)
                    {
                      const double v = hIn->GetBinContent(ix, iy) / sum;
                      h->SetBinContent(ix, iy, v);
                      h->SetBinError  (ix, iy, 0.0);
                    }
                  }
      
                  return h;
                };
      
                auto normalizeRecoFeedIn = [&](TH2* hIn, const char* newName)->TH2*
                {
                  if (!hIn) return nullptr;
      
                  TH2* h = CloneTH2(hIn, newName);
                  if (!h) return nullptr;
      
                  h->SetDirectory(nullptr);
                  EnsureSumw2(h);
                  h->Reset("ICES");
      
                  for (int ix = 1; ix <= hIn->GetNbinsX(); ++ix)
                  {
                    double sum = 0.0;
                    for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy) sum += hIn->GetBinContent(ix, iy);
                    if (!(sum > 0.0)) continue;
      
                    for (int iy = 1; iy <= hIn->GetNbinsY(); ++iy)
                    {
                      const double v = hIn->GetBinContent(ix, iy) / sum;
                      h->SetBinContent(ix, iy, v);
                      h->SetBinError  (ix, iy, 0.0);
                    }
                  }
      
                  return h;
                };
      
                auto buildResponseMeanRecoOverTruth = [&](TH2* hResp, const char* newName)->TH1D*
                {
                  if (!hResp) return nullptr;
      
                  TH1D* h = makeHistFromAxis(hResp->GetXaxis(), newName, "Mean(reco / truth)");
                  if (!h) return nullptr;
      
                  for (int ix = 1; ix <= hResp->GetNbinsX(); ++ix)
                  {
                    const double xTruth = hResp->GetXaxis()->GetBinCenter(ix);
                    if (!(xTruth > 0.0)) continue;
      
                    double sumW = 0.0;
                    double sumR = 0.0;
                    double sumR2 = 0.0;
      
                    for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy)
                    {
                      const double w = hResp->GetBinContent(ix, iy);
                      if (!(w > 0.0)) continue;
      
                      const double xReco = hResp->GetYaxis()->GetBinCenter(iy);
                      const double r = xReco / xTruth;
      
                      sumW  += w;
                      sumR  += w * r;
                      sumR2 += w * r * r;
                    }
      
                    if (!(sumW > 0.0)) continue;
      
                    const double mean = sumR / sumW;
                    const double var  = std::max(0.0, sumR2 / sumW - mean * mean);
                    const double err  = std::sqrt(var / sumW);
      
                    h->SetBinContent(ix, mean);
                    h->SetBinError  (ix, err);
                  }
      
                  return h;
                };
      
                auto buildResponseWidthRecoOverTruth = [&](TH2* hResp, const char* newName)->TH1D*
                {
                  if (!hResp) return nullptr;
      
                  TH1D* h = makeHistFromAxis(hResp->GetXaxis(), newName, "RMS(reco / truth)");
                  if (!h) return nullptr;
      
                  for (int ix = 1; ix <= hResp->GetNbinsX(); ++ix)
                  {
                    const double xTruth = hResp->GetXaxis()->GetBinCenter(ix);
                    if (!(xTruth > 0.0)) continue;
      
                    double sumW = 0.0;
                    double sumR = 0.0;
                    double sumR2 = 0.0;
      
                    for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy)
                    {
                      const double w = hResp->GetBinContent(ix, iy);
                      if (!(w > 0.0)) continue;
      
                      const double xReco = hResp->GetYaxis()->GetBinCenter(iy);
                      const double r = xReco / xTruth;
      
                      sumW  += w;
                      sumR  += w * r;
                      sumR2 += w * r * r;
                    }
      
                    if (!(sumW > 0.0)) continue;
      
                    const double mean = sumR / sumW;
                    const double var  = std::max(0.0, sumR2 / sumW - mean * mean);
                    const double rms  = std::sqrt(var);
                    const double err  = (sumW > 1.0) ? (rms / std::sqrt(2.0 * (sumW - 1.0))) : 0.0;
      
                    h->SetBinContent(ix, rms);
                    h->SetBinError  (ix, err);
                  }
      
                  return h;
                };
      
                auto buildResponseDiagonalFraction = [&](TH2* hResp, const char* newName)->TH1D*
                {
                  if (!hResp) return nullptr;
      
                  TH1D* h = makeHistFromAxis(hResp->GetXaxis(), newName, "Same-bin fraction");
                  if (!h) return nullptr;
      
                  for (int ix = 1; ix <= hResp->GetNbinsX(); ++ix)
                  {
                    const double xLo = hResp->GetXaxis()->GetBinLowEdge(ix);
                    const double xHi = hResp->GetXaxis()->GetBinUpEdge(ix);
      
                    int iyMatch = -1;
                    for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy)
                    {
                      const double yLo = hResp->GetYaxis()->GetBinLowEdge(iy);
                      const double yHi = hResp->GetYaxis()->GetBinUpEdge(iy);
                      if (std::fabs(yLo - xLo) < 1e-9 && std::fabs(yHi - xHi) < 1e-9)
                      {
                        iyMatch = iy;
                        break;
                      }
                    }
      
                    double sum = 0.0;
                    for (int iy = 1; iy <= hResp->GetNbinsY(); ++iy) sum += hResp->GetBinContent(ix, iy);
                    if (!(sum > 0.0) || iyMatch < 0) continue;
      
                    const double frac = hResp->GetBinContent(ix, iyMatch) / sum;
                    const double err  = std::sqrt(std::max(0.0, frac * (1.0 - frac) / sum));
      
                    h->SetBinContent(ix, frac);
                    h->SetBinError  (ix, err);
                  }
      
                  return h;
                };
      
                auto findExactBinByEdges = [&](const TAxis* ax, double lo, double hi)->int
                {
                  if (!ax) return -1;
                  for (int ib = 1; ib <= ax->GetNbins(); ++ib)
                  {
                    const double bLo = ax->GetBinLowEdge(ib);
                    const double bHi = ax->GetBinUpEdge(ib);
                    if (std::fabs(bLo - lo) < 1e-9 && std::fabs(bHi - hi) < 1e-9) return ib;
                  }
                  return -1;
                };
      
                auto styleHist = [&](TH1* h, int color, int marker)->void
                {
                  if (!h) return;
                  h->SetLineColor(color);
                  h->SetMarkerColor(color);
                  h->SetMarkerStyle(marker);
                  h->SetMarkerSize(0.95);
                  h->SetLineWidth(2);
                  h->SetTitle("");
                };
      
                auto prepPad = [&]()->void
                {
                  if (!gPad) return;
                  gPad->SetLeftMargin(0.12);
                  gPad->SetRightMargin(0.05);
                  gPad->SetBottomMargin(0.12);
                  gPad->SetTopMargin(0.10);
                  gPad->SetTicks(1, 1);
                };
      
                auto drawPadTitle = [&](const string& txt)->void
                {
                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(22);
                  tx.SetTextSize(0.050);
                  tx.DrawLatex(0.50, 0.95, txt.c_str());
                };
      
                auto drawCanvasTitle = [&](const string& txt, double size)->void
                {
                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(22);
                  tx.SetTextSize(size);
                  tx.DrawLatex(0.50, 0.985, txt.c_str());
                };
      
                auto drawCanvasNote = [&](const vector<string>& lines, double x, double y, int align, double size)->void
                {
                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(align);
                  tx.SetTextSize(size);
                  double yy = y;
                  for (const auto& line : lines)
                  {
                    tx.DrawLatex(x, yy, line.c_str());
                    yy -= 0.045;
                  }
                };
      
                auto drawMissingPad = [&](const string& yTitle, const string& msg)->void
                {
                  TH1F frame("frame_missing_phoQA", "", 1, 10.0, 35.0);
                  frame.SetMinimum(0.0);
                  frame.SetMaximum(1.0);
                  frame.SetTitle("");
                  frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                  frame.GetYaxis()->SetTitle(yTitle.c_str());
                  frame.Draw("axis");
      
                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(22);
                  tx.SetTextSize(0.050);
                  tx.DrawLatex(0.50, 0.50, msg.c_str());
                };
      
                auto setRatioRange = [&](TH1* h, double yMinDefault, double yMaxDefault)->void
                {
                  if (!h) return;
      
                  double yMin =  1e99;
                  double yMax = -1e99;
                  bool have = false;
      
                  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                  {
                    const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                    const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                    if (lo < 10.0 || hi > 35.0) continue;
      
                    const double y  = h->GetBinContent(ib);
                    const double ey = h->GetBinError(ib);
                    if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                    if (y <= 0.0 && ey <= 0.0) continue;
      
                    have = true;
                    yMin = std::min(yMin, y - ey);
                    yMax = std::max(yMax, y + ey);
                  }
      
                  if (!have || !(yMin < yMax))
                  {
                    h->GetYaxis()->SetRangeUser(yMinDefault, yMaxDefault);
                    return;
                  }
      
                  const double pad = std::max(0.08, 0.18 * (yMax - yMin));
                  yMin = std::min(yMin, 1.0) - pad;
                  yMax = std::max(yMax, 1.0) + pad;
                  if (yMin < 0.0) yMin = 0.0;
                  if (yMax <= yMin) yMax = yMin + 0.5;
                  h->GetYaxis()->SetRangeUser(yMin, yMax);
                };
      
                const bool buildPhotonQA = gApplyPurityCorrectionForUnfolding;
                const string phoQADir       = JoinPath(phoDir, "QA");
                const string qaMasterDir    = JoinPath(phoQADir, "01_master");
                const string qaBudgetDir    = JoinPath(phoQADir, "02_budget");
                const string qaResponseDir  = JoinPath(phoQADir, "03_response");
                const string qaSupportDir   = JoinPath(phoQADir, "04_support");
      
                if (buildPhotonQA)
                {
                  EnsureDir(phoQADir);
                  EnsureDir(qaMasterDir);
                  EnsureDir(qaBudgetDir);
                  EnsureDir(qaResponseDir);
                  EnsureDir(qaSupportDir);
                }
      
                if (hPhoResp_measXtruth)
                {
                  TCanvas c("c_pho_respT","c_pho_respT", 900, 750);
                  ApplyCanvasMargins2D(c);
                  c.SetLogz();
      
                  hPhoResp_measXtruth->SetTitle("");
                  hPhoResp_measXtruth->GetXaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                  hPhoResp_measXtruth->GetYaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
      
                  const auto& recoEdges  = kUnfoldRecoPtEdges;
                  const auto& truthEdges = kUnfoldTruthPtEdges;
      
                  const double recoMin     = (!recoEdges.empty() ? recoEdges.front() : 8.0);
                  const double recoUFHi    = (recoEdges.size() >= 2 ? recoEdges[1] : 10.0);
                  const double recoOFLo    = (recoEdges.size() >= 2 ? recoEdges[recoEdges.size() - 2] : 35.0);
                  const double recoDrawMax = recoOFLo;
      
                  const double truthMin   = (!truthEdges.empty() ? truthEdges.front() : 5.0);
                  const double truthUFMid = (truthEdges.size() >= 2 ? truthEdges[1] : 8.0);
                  const double truthUFHi  = (truthEdges.size() >= 3 ? truthEdges[2] : 10.0);
                  const double truthOFLo  = (truthEdges.size() >= 2 ? truthEdges[truthEdges.size() - 2] : 35.0);
                  const double truthMax   = (!truthEdges.empty() ? truthEdges.back() : 40.0);
      
                  hPhoResp_measXtruth->GetXaxis()->SetRangeUser(recoMin, recoDrawMax);
                  hPhoResp_measXtruth->GetYaxis()->SetRangeUser(truthMin, truthMax);
                  hPhoResp_measXtruth->Draw("colz");
                  if (gPad) gPad->Update();
      
                  TLine lRecoUF(recoUFHi, truthMin, recoUFHi, truthMax);
                  TLine lRecoOF(recoOFLo, truthMin, recoOFLo, truthMax);
                  TLine lTruthUF1(recoMin, truthUFMid, recoDrawMax, truthUFMid);
                  TLine lTruthUF2(recoMin, truthUFHi, recoDrawMax, truthUFHi);
                  TLine lTruthOF(recoMin, truthOFLo, recoDrawMax, truthOFLo);
      
                  lRecoUF.SetLineColor(kBlack);
                  lRecoOF.SetLineColor(kBlack);
                  lTruthUF1.SetLineColor(kBlack);
                  lTruthUF2.SetLineColor(kBlack);
                  lTruthOF.SetLineColor(kBlack);
      
                  lRecoUF.SetLineStyle(2);
                  lRecoOF.SetLineStyle(2);
                  lTruthUF1.SetLineStyle(2);
                  lTruthUF2.SetLineStyle(2);
                  lTruthOF.SetLineStyle(2);
      
                  lRecoUF.SetLineWidth(2);
                  lRecoOF.SetLineWidth(2);
                  lTruthUF1.SetLineWidth(2);
                  lTruthUF2.SetLineWidth(2);
                  lTruthOF.SetLineWidth(2);
      
                  if (recoUFHi > recoMin) lRecoUF.Draw("same");
                  if (recoOFLo > recoMin) lRecoOF.Draw("same");
                  if (truthUFMid > truthMin) lTruthUF1.Draw("same");
                  if (truthUFHi > truthMin) lTruthUF2.Draw("same");
                  if (truthOFLo > truthMin) lTruthOF.Draw("same");
      
                  DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsSim), 0.034, 0.045);
                  DrawLatexLines(0.14,0.88, { "SIM photon response (transpose)", "reco #rightarrow truth axis order" }, 0.030, 0.040);
                  DrawLatexLines(0.52,0.36,
                                 {
                                   "UF / OF support bins:",
                                   "reco: 8-10 = UF, 35-40 = OF",
                                   "truth: 5-8 and 8-10 = UF, 35-40 = OF"
                                 },
                                 0.026, 0.035);
      
                  c.Modified();
                  c.Update();
                  SaveCanvas(c, JoinPath(phoDir, "pho_response_recoVsTruth.png"));
                }
      
                if (hPhoUnfoldTruth)
                {
                  TH1* hRecoShape = CloneTH1(hPhoRecoData, "hPhoRecoData_forOverlayLog");
                  TH1* hUnfShape  = CloneTH1(hPhoUnfoldTruth, "hPhoUnfoldTruth_forOverlayLog");
      
                  if (hRecoShape && hUnfShape)
                  {
                    hRecoShape->SetDirectory(nullptr);
                    hUnfShape->SetDirectory(nullptr);
                    EnsureSumw2(hRecoShape);
                    EnsureSumw2(hUnfShape);
      
                    hRecoShape->SetLineColor(2);
                    hRecoShape->SetMarkerColor(2);
                    hRecoShape->SetMarkerStyle(20);
                    hRecoShape->SetLineWidth(2);
      
                    hUnfShape->SetLineColor(1);
                    hUnfShape->SetMarkerColor(1);
                    hUnfShape->SetMarkerStyle(24);
                    hUnfShape->SetLineWidth(2);
      
                    const double xPlotMin = 10.0;
                    const double xPlotMax = 35.0;
      
                    auto scanMaxInRange = [&](TH1* h)->double
                    {
                      double maxY = 0.0;
                      if (!h) return maxY;
      
                      const int nb = h->GetNbinsX();
                      for (int ib = 1; ib <= nb; ++ib)
                      {
                        const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                        const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                        if (lo < xPlotMin || hi > xPlotMax) continue;
      
                        const double y  = h->GetBinContent(ib);
                        const double ey = h->GetBinError(ib);
                        if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                        if (y <= 0.0 && ey <= 0.0) continue;
      
                        const double v = y + ey;
                        if (v > maxY) maxY = v;
                      }
                      return maxY;
                    };
      
                    auto scanMinPosInRange = [&](TH1* h)->double
                    {
                      double minY = 1e99;
                      if (!h) return minY;
      
                      const int nb = h->GetNbinsX();
                      for (int ib = 1; ib <= nb; ++ib)
                      {
                        const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                        const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                        if (lo < xPlotMin || hi > xPlotMax) continue;
      
                        const double y = h->GetBinContent(ib);
                        if (!std::isfinite(y)) continue;
                        if (y > 0.0 && y < minY) minY = y;
                      }
                      return minY;
                    };
      
                    const double maxvTop = std::max(scanMaxInRange(hRecoShape), scanMaxInRange(hUnfShape));
                    double minPos = std::min(scanMinPosInRange(hRecoShape), scanMinPosInRange(hUnfShape));
                    if (!(minPos < 1e98)) minPos = 1e-3;
      
                    TH1* hRatio = makeRatioOnRef(hPhoUnfoldTruth, hPhoRecoData, hPhoRecoData, "hPho_afterOverBefore_ratio");
                    if (hRatio)
                    {
                      hRatio->SetDirectory(nullptr);
                      EnsureSumw2(hRatio);
                    }
      
                    TCanvas c("c_pho_unf_logy","c_pho_unf_logy", 900, 700);
                    TPad* pTop = new TPad("pTop_pho_unf_logy", "pTop_pho_unf_logy", 0.0, 0.30, 1.0, 1.0);
                    TPad* pBot = new TPad("pBot_pho_unf_logy", "pBot_pho_unf_logy", 0.0, 0.00, 1.0, 0.30);
      
                    pTop->SetLeftMargin(0.12);
                    pTop->SetRightMargin(0.04);
                    pTop->SetTopMargin(0.08);
                    pTop->SetBottomMargin(0.02);
                    pTop->SetTicks(1, 1);
                    pTop->SetLogy(1);
      
                    pBot->SetLeftMargin(0.12);
                    pBot->SetRightMargin(0.04);
                    pBot->SetTopMargin(0.02);
                    pBot->SetBottomMargin(0.30);
                    pBot->SetTicks(1, 1);
                    pBot->SetGridy(1);
      
                    pTop->Draw();
                    pBot->Draw();
      
                    pTop->cd();
                    hRecoShape->SetTitle("");
                    hRecoShape->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                    hUnfShape->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                    hRecoShape->GetXaxis()->SetTitle("");
                    hRecoShape->GetXaxis()->SetLabelSize(0.0);
                    hRecoShape->GetXaxis()->SetTitleSize(0.0);
                    hRecoShape->GetYaxis()->SetTitle("Counts");
                    hRecoShape->GetYaxis()->SetTitleSize(0.055);
                    hRecoShape->GetYaxis()->SetTitleOffset(0.95);
                    hRecoShape->GetYaxis()->SetLabelSize(0.045);
                    hRecoShape->SetMinimum(std::max(minPos * 0.5, 1e-6));
                    hRecoShape->SetMaximum((maxvTop > 0.0) ? (maxvTop * 3.0) : 1.0);
                    hRecoShape->Draw("E1");
                    hUnfShape->Draw("E1 same");
      
                    TLegend leg(0.55,0.76,0.92,0.90);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.032);
                    leg.AddEntry(hRecoShape, "PP DATA (reco)", "lep");
                    leg.AddEntry(hUnfShape,  TString::Format("Unfolded (truth), Bayes it=%d", kBayesIterPho).Data(), "lep");
                    leg.Draw();
      
                    DrawLatexLines(0.14,0.92,
                                   { "Photon 1D Unfolding, N_{#gamma}(p_{T}^{#gamma}), Photon 4 GeV + MBD NS #geq 1" },
                                   0.034, 0.045);
      
                    pBot->cd();
                    if (hRatio)
                    {
                      hRatio->SetTitle("");
                      hRatio->SetMarkerStyle(20);
                      hRatio->SetMarkerSize(0.95);
                      hRatio->SetLineWidth(2);
                      hRatio->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                      hRatio->GetYaxis()->SetTitle("After / before unfolding");
                      hRatio->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                      setRatioRange(hRatio, 0.8, 1.2);
                      hRatio->GetXaxis()->SetTitleSize(0.12);
                      hRatio->GetXaxis()->SetLabelSize(0.11);
                      hRatio->GetXaxis()->SetTitleOffset(1.00);
                      hRatio->GetYaxis()->SetTitleSize(0.10);
                      hRatio->GetYaxis()->SetLabelSize(0.09);
                      hRatio->GetYaxis()->SetTitleOffset(0.55);
                      hRatio->GetYaxis()->SetNdivisions(505);
                      hRatio->Draw("E1");
      
                      TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");
                    }
                    else
                    {
                      TH1F frame("frame_pho_ratio","", 1, xPlotMin, xPlotMax);
                      frame.SetMinimum(0.8);
                      frame.SetMaximum(1.2);
                      frame.SetTitle("");
                      frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                      frame.GetYaxis()->SetTitle("After / before unfolding");
                      frame.GetXaxis()->SetTitleSize(0.12);
                      frame.GetXaxis()->SetLabelSize(0.11);
                      frame.GetXaxis()->SetTitleOffset(1.00);
                      frame.GetYaxis()->SetTitleSize(0.10);
                      frame.GetYaxis()->SetLabelSize(0.09);
                      frame.GetYaxis()->SetTitleOffset(0.55);
                      frame.GetYaxis()->SetNdivisions(505);
                      frame.Draw("axis");
      
                      TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");
                    }
      
                    c.cd();
                    c.Modified();
                    c.Update();
                    SaveCanvas(c, JoinPath(phoDir, "pho_unfolded_truth_pTgamma_overlay_logy.png"));
      
                    if (hRatio) delete hRatio;
                    delete pTop;
                    delete pBot;
                  }
      
                  if (hRecoShape) delete hRecoShape;
                  if (hUnfShape) delete hUnfShape;
                }
      
                if (buildPhotonQA)
                {
                  TH1* hRecoDataQA   = cloneOrMapToRef(hPhoRecoData, hPhoRecoData, "hPhoRecoData_QA");
                  TH1* hRecoSimQA    = cloneOrMapToRef(hPhoRecoSim,  hPhoRecoData, "hPhoRecoSim_QA");
                  TH1* hTruthSimQA   = cloneOrMapToRef(hPhoTruthSim, hPhoRecoData, "hPhoTruthSim_onReco_QA");
                  TH1* hMissesSimQA  = cloneOrMapToRef(hPhoTruthMissesSim_in, hPhoRecoData, "hPhoTruthMisses_onReco_QA");
                  TH1* hUnfoldDataQA = cloneOrMapToRef(hPhoUnfoldTruth, hPhoRecoData, "hPhoUnfoldTruth_onReco_QA");
      
                  TH1* hTruthOverRecoSim = makeRatioOnRef(hPhoTruthSim, hPhoRecoSim, hPhoRecoData, "hPho_truthOverRecoSim_QA");
                  TH1* hRecoOverData     = makeRatioOnRef(hPhoRecoSim, hPhoRecoData, hPhoRecoData, "hPho_recoSimOverData_QA");
                  TH1* hTruthOverData    = makeRatioOnRef(hPhoTruthSim, hPhoRecoData, hPhoRecoData, "hPho_truthSimOverData_QA");
      
                  TH1* hFakeFracSim = makeRatioOnRef(hPhoRecoFakesSim_in, hPhoRecoSim, hPhoRecoData, "hPho_fakeFractionSim_QA");
                  TH1* hPuritySim   = makeOneMinus(hFakeFracSim, "hPho_puritySim_QA");
                  TH1* hMissFracSim = makeRatioOnRef(hPhoTruthMissesSim_in, hPhoTruthSim, hPhoRecoData, "hPho_missFractionSim_QA");
                  TH1* hEffSim      = makeOneMinus(hMissFracSim, "hPho_efficiencySim_QA");
                  TH1* hAfterOverBefore = makeRatioOnRef(hPhoUnfoldTruth, hPhoRecoData, hPhoRecoData, "hPho_afterOverBefore_QA");
      
                  styleHist(hRecoDataQA,   kBlack,    20);
                  styleHist(hRecoSimQA,    kBlue + 1, 24);
                  styleHist(hTruthSimQA,   kRed + 1,  25);
                  styleHist(hTruthOverRecoSim, kBlack,    20);
                  styleHist(hRecoOverData,     kBlue + 1, 20);
                  styleHist(hTruthOverData,    kRed + 1,  20);
                  styleHist(hFakeFracSim,      kRed + 1,  20);
                  styleHist(hPuritySim,        kBlue + 1, 20);
                  styleHist(hMissFracSim,      kRed + 1,  20);
                  styleHist(hEffSim,           kBlue + 1, 20);
                  styleHist(hAfterOverBefore,  kBlack,    20);
      
                  if (hRecoDataQA && hRecoSimQA && hTruthSimQA && hTruthOverRecoSim && hRecoOverData && hTruthOverData)
                  {
                    TCanvas c("c_pho_masterQA", "c_pho_masterQA", 1500, 1050);
                    c.Divide(2, 2, 0.001, 0.001);
      
                    c.cd(1);
                    prepPad();
                    if (gPad) gPad->SetLogy(1);
      
                    double maxTop = 0.0;
                    auto scanMasterMax = [&](TH1* h)->void
                    {
                      if (!h) return;
                      for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                      {
                        const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                        const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                        if (lo < 10.0 || hi > 35.0) continue;
                        const double v = h->GetBinContent(ib) + h->GetBinError(ib);
                        if (std::isfinite(v) && v > maxTop) maxTop = v;
                      }
                    };
                    scanMasterMax(hRecoDataQA);
                    scanMasterMax(hRecoSimQA);
                    scanMasterMax(hTruthSimQA);
      
                    double minPos = 1e99;
                    auto scanMasterMinPos = [&](TH1* h)->void
                    {
                      if (!h) return;
                      for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                      {
                        const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                        const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                        if (lo < 10.0 || hi > 35.0) continue;
                        const double v = h->GetBinContent(ib);
                        if (std::isfinite(v) && v > 0.0 && v < minPos) minPos = v;
                      }
                    };
                    scanMasterMinPos(hRecoDataQA);
                    scanMasterMinPos(hRecoSimQA);
                    scanMasterMinPos(hTruthSimQA);
                    if (!(minPos < 1e98)) minPos = 1e-3;
      
                    hRecoDataQA->SetMinimum(std::max(0.5 * minPos, 1e-6));
                    hRecoDataQA->SetMaximum((maxTop > 0.0) ? (3.0 * maxTop) : 1.0);
                    hRecoDataQA->GetXaxis()->SetRangeUser(10.0, 35.0);
                    hRecoDataQA->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    hRecoDataQA->GetYaxis()->SetTitle("Counts");
                    hRecoDataQA->Draw("E1");
                    hRecoSimQA->GetXaxis()->SetRangeUser(10.0, 35.0);
                    hRecoSimQA->Draw("E1 same");
                    hTruthSimQA->GetXaxis()->SetRangeUser(10.0, 35.0);
                    hTruthSimQA->Draw("E1 same");
                    drawPadTitle("DATA reco vs SIM reco vs SIM truth");
      
                    TLegend leg(0.52, 0.70, 0.90, 0.88);
                    leg.SetBorderSize(0);
                    leg.SetFillStyle(0);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.040);
                    leg.AddEntry(hRecoDataQA, "DATA reco", "pe");
                    leg.AddEntry(hRecoSimQA,  "SIM reco",  "pe");
                    leg.AddEntry(hTruthSimQA, "SIM truth (bin-mapped)", "pe");
                    leg.Draw();
      
                    c.cd(2);
                    prepPad();
                    hTruthOverRecoSim->GetXaxis()->SetRangeUser(10.0, 35.0);
                    hTruthOverRecoSim->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    hTruthOverRecoSim->GetYaxis()->SetTitle("SIM truth / SIM reco");
                    setRatioRange(hTruthOverRecoSim, 0.5, 2.5);
                    hTruthOverRecoSim->Draw("E1");
                    {
                      TLine l1(10.0, 1.0, 35.0, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");
                    }
                    drawPadTitle("SIM lever-arm");
      
                    c.cd(3);
                    prepPad();
                    hRecoOverData->GetXaxis()->SetRangeUser(10.0, 35.0);
                    hRecoOverData->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    hRecoOverData->GetYaxis()->SetTitle("SIM reco / DATA reco");
                    setRatioRange(hRecoOverData, 0.5, 1.5);
                    hRecoOverData->Draw("E1");
                    {
                      TLine l1(10.0, 1.0, 35.0, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");
                    }
                    drawPadTitle("Reco data/MC check");
      
                    c.cd(4);
                    prepPad();
                    hTruthOverData->GetXaxis()->SetRangeUser(10.0, 35.0);
                    hTruthOverData->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    hTruthOverData->GetYaxis()->SetTitle("SIM truth / DATA reco");
                    setRatioRange(hTruthOverData, 0.5, 2.5);
                    hTruthOverData->Draw("E1");
                    {
                      TLine l1(10.0, 1.0, 35.0, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");
                    }
                    drawPadTitle("Overall truth/data gap");
      
                    c.cd();
                    drawCanvasTitle("Photon master QA: yields and key ratios", 0.028);
                    drawCanvasNote(
                      {
                        "Displayed bins: 10-12, 12-14, 14-16, 16-18, 18-20, 20-22, 22-24, 24-26, 26-35 GeV",
                        "SIM truth is mapped onto the reco analysis axis when needed"
                      },
                      0.97, 0.05, 31, 0.020
                    );
                    SaveCanvas(c, JoinPath(qaMasterDir, "pho_master_yields_and_ratios.png"));
                  }
      
                  if (hFakeFracSim && hPuritySim && hMissFracSim && hEffSim && hTruthOverRecoSim && hAfterOverBefore)
                  {
                    TCanvas c("c_pho_budgetQA", "c_pho_budgetQA", 1700, 1050);
                    c.Divide(3, 2, 0.001, 0.001);
      
                    vector<pair<TH1*, string>> pads =
                    {
                      { hFakeFracSim,     "Fake fraction" },
                      { hPuritySim,       "Purity = 1 - fake/reco" },
                      { hMissFracSim,     "Miss fraction" },
                      { hEffSim,          "Efficiency = 1 - miss/truth" },
                      { hTruthOverRecoSim,"SIM truth / SIM reco" },
                      { hAfterOverBefore, "DATA after / before unfolding" }
                    };
      
                    for (int ipad = 0; ipad < 6; ++ipad)
                    {
                      c.cd(ipad + 1);
                      prepPad();
      
                      TH1* h = pads[(std::size_t)ipad].first;
                      if (!h)
                      {
                        drawMissingPad("Value", "MISSING");
                        continue;
                      }
      
                      h->GetXaxis()->SetRangeUser(10.0, 35.0);
                      h->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                      h->GetYaxis()->SetTitle("Value");
      
                      if (ipad < 4)
                      {
                        h->GetYaxis()->SetRangeUser(0.0, 1.2);
                      }
                      else
                      {
                        setRatioRange(h, 0.5, 2.5);
                      }
      
                      h->Draw("E1");
      
                      if (ipad != 0 && ipad != 2)
                      {
                        TLine l1(10.0, 1.0, 35.0, 1.0);
                        l1.SetLineStyle(2);
                        l1.SetLineWidth(2);
                        l1.Draw("same");
                      }
      
                      drawPadTitle(pads[(std::size_t)ipad].second);
                    }
      
                    c.cd();
                    drawCanvasTitle("Photon correction-budget QA", 0.028);
                    drawCanvasNote(
                      {
                        "Top row: reco fake contamination and truth-matching loss terms",
                        "Bottom-right compares the observed data correction to the SIM truth/reco lever-arm"
                      },
                      0.97, 0.05, 31, 0.020
                    );
                    SaveCanvas(c, JoinPath(qaBudgetDir, "pho_correction_budget_summary.png"));
                  }
      
                  if (hPhoRespSim)
                  {
                    TH2* hRespTruthNorm = normalizeResponseByTruth(hPhoRespSim, "hPhoResp_truthNorm_QA");
                    TH2* hRespRecoNorm  = normalizeResponseByReco (hPhoRespSim, "hPhoResp_recoNorm_QA");
                    TH1D* hRespMean     = buildResponseMeanRecoOverTruth(hPhoRespSim, "hPhoResp_meanRecoOverTruth_QA");
                    TH1D* hRespWidth    = buildResponseWidthRecoOverTruth(hPhoRespSim, "hPhoResp_widthRecoOverTruth_QA");
                    TH1D* hRespDiagFrac = buildResponseDiagonalFraction(hPhoRespSim, "hPhoResp_diagFrac_QA");
      
                    styleHist(hRespMean,     kBlack,    20);
                    styleHist(hRespWidth,    kBlue + 1, 20);
                    styleHist(hRespDiagFrac, kRed + 1,  20);
      
                    TCanvas c("c_pho_responseShapeQA", "c_pho_responseShapeQA", 1750, 1050);
                    c.Divide(3, 2, 0.001, 0.001);
      
                    c.cd(1);
                    prepPad();
                    if (gPad) gPad->SetLogz(1);
                    TH2* hRawDraw = CloneTH2(hPhoRespSim, "hPhoResp_rawDraw_QA");
                    if (hRawDraw)
                    {
                      hRawDraw->SetDirectory(nullptr);
                      hRawDraw->SetTitle("");
                      hRawDraw->GetXaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
                      hRawDraw->GetYaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                      double minPos = 1e99;
                      for (int ix = 1; ix <= hRawDraw->GetNbinsX(); ++ix)
                      {
                        for (int iy = 1; iy <= hRawDraw->GetNbinsY(); ++iy)
                        {
                          const double z = hRawDraw->GetBinContent(ix, iy);
                          if (std::isfinite(z) && z > 0.0 && z < minPos) minPos = z;
                        }
                      }
                      if (minPos < 1e98) hRawDraw->SetMinimum(std::max(0.5 * minPos, 1e-6));
                      hRawDraw->Draw("colz");
                      drawPadTitle("Raw response counts");
                    }
      
                    c.cd(2);
                    prepPad();
                    if (hRespTruthNorm)
                    {
                      hRespTruthNorm->SetTitle("");
                      hRespTruthNorm->GetXaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
                      hRespTruthNorm->GetYaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                      hRespTruthNorm->SetMinimum(0.0);
                      hRespTruthNorm->SetMaximum(1.0);
                      hRespTruthNorm->Draw("colz");
                      drawPadTitle("Truth-bin normalized: P(reco|truth)");
                    }
                    else
                    {
                      drawMissingPad("Probability", "MISSING");
                    }
      
                    c.cd(3);
                    prepPad();
                    if (hRespRecoNorm)
                    {
                      hRespRecoNorm->SetTitle("");
                      hRespRecoNorm->GetXaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
                      hRespRecoNorm->GetYaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                      hRespRecoNorm->SetMinimum(0.0);
                      hRespRecoNorm->SetMaximum(1.0);
                      hRespRecoNorm->Draw("colz");
                      drawPadTitle("Reco-bin normalized: P(truth|reco)");
                    }
                    else
                    {
                      drawMissingPad("Probability", "MISSING");
                    }
      
                    c.cd(4);
                    prepPad();
                    if (hRespMean)
                    {
                      hRespMean->GetXaxis()->SetRangeUser(5.0, 40.0);
                      hRespMean->GetYaxis()->SetRangeUser(0.6, 1.4);
                      hRespMean->Draw("E1");
                      TLine l1(5.0, 1.0, 40.0, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");
                      drawPadTitle("Mean reco/truth vs truth p_{T}");
                    }
                    else
                    {
                      drawMissingPad("Mean(reco / truth)", "MISSING");
                    }
      
                    c.cd(5);
                    prepPad();
                    if (hRespWidth)
                    {
                      hRespWidth->GetXaxis()->SetRangeUser(5.0, 40.0);
                      hRespWidth->GetYaxis()->SetRangeUser(0.0, 0.6);
                      hRespWidth->Draw("E1");
                      drawPadTitle("RMS(reco/truth) vs truth p_{T}");
                    }
                    else
                    {
                      drawMissingPad("RMS(reco / truth)", "MISSING");
                    }
      
                    c.cd(6);
                    prepPad();
                    if (hRespDiagFrac)
                    {
                      hRespDiagFrac->GetXaxis()->SetRangeUser(5.0, 40.0);
                      hRespDiagFrac->GetYaxis()->SetRangeUser(0.0, 1.05);
                      hRespDiagFrac->Draw("E1");
                      drawPadTitle("Same-bin fraction vs truth p_{T}");
                    }
                    else
                    {
                      drawMissingPad("Same-bin fraction", "MISSING");
                    }
      
                    c.cd();
                    drawCanvasTitle("Photon response-shape QA", 0.028);
                    drawCanvasNote(
                      {
                        "Top row shows the raw response and its truth/reco normalizations",
                        "Bottom row compresses the matrix into migration bias, width, and diagonal fraction"
                      },
                      0.97, 0.05, 31, 0.020
                    );
                    SaveCanvas(c, JoinPath(qaResponseDir, "pho_responseShape_summary.png"));
      
                    if (hRawDraw) delete hRawDraw;
                    if (hRespTruthNorm) delete hRespTruthNorm;
                    if (hRespRecoNorm) delete hRespRecoNorm;
                    if (hRespMean) delete hRespMean;
                    if (hRespWidth) delete hRespWidth;
                    if (hRespDiagFrac) delete hRespDiagFrac;
                  }
      
                  if (hPhoResp_measXtruth)
                  {
                    TH2* hFeedIn = normalizeRecoFeedIn(hPhoResp_measXtruth, "hPhoResp_recoFeedIn_QA");
                    TH1* hFeed58 = CloneTH1(hPhoRecoData, "hPho_feedFromTruth5to8_QA");
                    TH1* hFeed810 = CloneTH1(hPhoRecoData, "hPho_feedFromTruth8to10_QA");
                    TH1* hFeedSupport = CloneTH1(hPhoRecoData, "hPho_feedFromTruthSupport_QA");
      
                    if (hFeed58)     { hFeed58->SetDirectory(nullptr);     EnsureSumw2(hFeed58);     hFeed58->Reset("ICES"); }
                    if (hFeed810)    { hFeed810->SetDirectory(nullptr);    EnsureSumw2(hFeed810);    hFeed810->Reset("ICES"); }
                    if (hFeedSupport){ hFeedSupport->SetDirectory(nullptr);EnsureSumw2(hFeedSupport);hFeedSupport->Reset("ICES"); }
      
                    const int iy58  = (hFeedIn ? findExactBinByEdges(hFeedIn->GetYaxis(), 5.0, 8.0)  : -1);
                    const int iy810 = (hFeedIn ? findExactBinByEdges(hFeedIn->GetYaxis(), 8.0, 10.0) : -1);
      
                    if (hFeedIn && hFeed58 && hFeed810 && hFeedSupport)
                    {
                      for (int ix = 1; ix <= hFeedIn->GetNbinsX(); ++ix)
                      {
                        const double recoLo = hFeedIn->GetXaxis()->GetBinLowEdge(ix);
                        const double recoHi = hFeedIn->GetXaxis()->GetBinUpEdge(ix);
                        const double recoCtr = hFeedIn->GetXaxis()->GetBinCenter(ix);
                        const int ibRef = hFeed58->GetXaxis()->FindBin(recoCtr);
                        if (ibRef < 1 || ibRef > hFeed58->GetNbinsX()) continue;
      
                        if (recoLo < 10.0 || recoHi > 35.0) continue;
      
                        double f58  = (iy58  > 0 ? hFeedIn->GetBinContent(ix, iy58)  : 0.0);
                        double f810 = (iy810 > 0 ? hFeedIn->GetBinContent(ix, iy810) : 0.0);
                        const double fSup = f58 + f810;
      
                        hFeed58->SetBinContent(ibRef, f58);
                        hFeed810->SetBinContent(ibRef, f810);
                        hFeedSupport->SetBinContent(ibRef, fSup);
                        hFeed58->SetBinError(ibRef, 0.0);
                        hFeed810->SetBinError(ibRef, 0.0);
                        hFeedSupport->SetBinError(ibRef, 0.0);
                      }
      
                      styleHist(hFeed58,      kBlue + 1, 24);
                      styleHist(hFeed810,     kRed + 1,  25);
                      styleHist(hFeedSupport, kBlack,    20);
      
                      TCanvas c("c_pho_supportQA", "c_pho_supportQA", 1600, 760);
                      c.Divide(2, 1, 0.001, 0.001);
      
                      c.cd(1);
                      prepPad();
                      if (hFeedIn)
                      {
                        hFeedIn->SetTitle("");
                        hFeedIn->GetXaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                        hFeedIn->GetYaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
                        hFeedIn->GetXaxis()->SetRangeUser(10.0, 35.0);
                        hFeedIn->GetYaxis()->SetRangeUser(5.0, 40.0);
                        hFeedIn->SetMinimum(0.0);
                        hFeedIn->SetMaximum(1.0);
                        hFeedIn->Draw("colz");
      
                        if (iy58 > 0)
                        {
                          TLine l(10.0, hFeedIn->GetYaxis()->GetBinUpEdge(iy58), 35.0, hFeedIn->GetYaxis()->GetBinUpEdge(iy58));
                          l.SetLineStyle(2);
                          l.SetLineWidth(2);
                          l.Draw("same");
                        }
                        if (iy810 > 0)
                        {
                          TLine l(10.0, hFeedIn->GetYaxis()->GetBinUpEdge(iy810), 35.0, hFeedIn->GetYaxis()->GetBinUpEdge(iy810));
                          l.SetLineStyle(2);
                          l.SetLineWidth(2);
                          l.Draw("same");
                        }
      
                        drawPadTitle("Reco-bin feed-in from truth bins");
                      }
                      else
                      {
                        drawMissingPad("Fraction", "MISSING");
                      }
      
                      c.cd(2);
                      prepPad();
                      if (hFeed58 && hFeed810 && hFeedSupport)
                      {
                        double maxY = 0.0;
                        auto scan = [&](TH1* h)->void
                        {
                          for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                          {
                            const double lo = h->GetXaxis()->GetBinLowEdge(ib);
                            const double hi = h->GetXaxis()->GetBinUpEdge(ib);
                            if (lo < 10.0 || hi > 35.0) continue;
                            const double v = h->GetBinContent(ib) + h->GetBinError(ib);
                            if (std::isfinite(v) && v > maxY) maxY = v;
                          }
                        };
                        scan(hFeed58);
                        scan(hFeed810);
                        scan(hFeedSupport);
      
                        hFeedSupport->GetXaxis()->SetRangeUser(10.0, 35.0);
                        hFeedSupport->GetXaxis()->SetTitle("p_{T}^{#gamma, reco} [GeV]");
                        hFeedSupport->GetYaxis()->SetTitle("Feed-in fraction");
                        hFeedSupport->GetYaxis()->SetRangeUser(0.0, (maxY > 0.0) ? (1.15 * maxY) : 0.20);
                        hFeedSupport->Draw("E1");
                        hFeed58->GetXaxis()->SetRangeUser(10.0, 35.0);
                        hFeed810->GetXaxis()->SetRangeUser(10.0, 35.0);
                        hFeed58->Draw("E1 same");
                        hFeed810->Draw("E1 same");
                        drawPadTitle("Low-p_{T}^{truth} support feed-in");
      
                        TLegend leg(0.48, 0.68, 0.90, 0.88);
                        leg.SetBorderSize(0);
                        leg.SetFillStyle(0);
                        leg.SetTextFont(42);
                        leg.SetTextSize(0.040);
                        leg.AddEntry(hFeedSupport, "Total support feed-in (5-10 GeV truth)", "pe");
                        leg.AddEntry(hFeed58,      "5-8 GeV truth", "pe");
                        leg.AddEntry(hFeed810,     "8-10 GeV truth", "pe");
                        leg.Draw();
                      }
                      else
                      {
                        drawMissingPad("Feed-in fraction", "MISSING");
                      }
      
                      c.cd();
                      drawCanvasTitle("Photon support-bin feed-in QA", 0.028);
                      drawCanvasNote(
                        {
                          "Left: each reco analysis bin normalized to unit truth-bin feed-in",
                          "Right: isolates how much of each displayed reco bin comes from the low-p_{T}^{truth} support bins"
                        },
                        0.97, 0.05, 31, 0.020
                      );
                      SaveCanvas(c, JoinPath(qaSupportDir, "pho_supportFeedin_summary.png"));
                    }
      
                    if (hFeedIn) delete hFeedIn;
                    if (hFeed58) delete hFeed58;
                    if (hFeed810) delete hFeed810;
                    if (hFeedSupport) delete hFeedSupport;
                  }
      
                  if (hAfterOverBefore) delete hAfterOverBefore;
                  if (hEffSim) delete hEffSim;
                  if (hMissFracSim) delete hMissFracSim;
                  if (hPuritySim) delete hPuritySim;
                  if (hFakeFracSim) delete hFakeFracSim;
                  if (hTruthOverData) delete hTruthOverData;
                  if (hRecoOverData) delete hRecoOverData;
                  if (hTruthOverRecoSim) delete hTruthOverRecoSim;
                  if (hUnfoldDataQA) delete hUnfoldDataQA;
                  if (hMissesSimQA) delete hMissesSimQA;
                  if (hTruthSimQA) delete hTruthSimQA;
                  if (hRecoSimQA) delete hRecoSimQA;
                  if (hRecoDataQA) delete hRecoDataQA;
                }
              }
      
              // Photon summary for unfolding RECO pT bins
              vector<string> phoSummary;
              phoSummary.push_back("Photon unfolding summary (PP DATA unfolded to truth pTgamma)");
              phoSummary.push_back(TString::Format("Method: RooUnfoldBayes, iterations=%d, error=kCovToy, Ntoys=%d", kBayesIterPho, kNToysPhoFinal).Data());
              phoSummary.push_back("");
      
              {
                  const auto& analysisRecoBins = UnfoldAnalysisRecoPtBins();
                  const int nPtUnf = (int)analysisRecoBins.size();
                  for (int i = 0; i < nPtUnf; ++i)
                  {
                    const PtBin& b = analysisRecoBins[i];
                    const double cen = 0.5 * (b.lo + b.hi);
      
                    const int ibTruth = (hPhoUnfoldTruth ? hPhoUnfoldTruth->GetXaxis()->FindBin(cen) : -1);
      
                    const string labCanon = TString::Format("%d-%d GeV", b.lo, b.hi).Data();
                    const string labTruth = (hPhoUnfoldTruth ? AxisBinLabel(hPhoUnfoldTruth->GetXaxis(), ibTruth, "GeV", 0) : "N/A");
      
                    const double N = (hPhoUnfoldTruth ? hPhoUnfoldTruth->GetBinContent(ibTruth) : 0.0);
                    const double E = (hPhoUnfoldTruth ? hPhoUnfoldTruth->GetBinError  (ibTruth) : 0.0);
      
                    phoSummary.push_back(
                      TString::Format("pT^gamma analysis=%s  -> truthBin=%s  N_gamma(unf)=%.6g ± %.6g",
                        labCanon.c_str(), labTruth.c_str(), N, E
                      ).Data()
                    );
                  }
                  phoSummary.push_back("");
                  phoSummary.push_back("NOTE: only the 9 analysis bins are summarized/plotted: 10-12, 12-14, 14-16, 16-18, 18-20, 20-22, 22-24, 24-26, 26-35 GeV.");
                  phoSummary.push_back("      The reco 8-10 GeV bin is retained only as an unfolding support underflow bin.");
                  phoSummary.push_back("      The truth 5-8 GeV and 8-10 GeV bins are retained only as truth underflow support bins.");
                  phoSummary.push_back("      The reco/truth 35-40 GeV bin is retained only as an overflow support bin.");
                  phoSummary.push_back("      RooUnfold still uses the full YAML unfolding axes; only the displayed analysis-bin outputs are restricted here.");
              }
      
            WriteTextFile(JoinPath(phoDir, "summary_photon_unfolding_bins.txt"), phoSummary);
      
            // Save photon ROOT outputs
            {
              const string outRoot = JoinPath(phoDir, "rooUnfold_photons.root");
              TFile f(outRoot.c_str(), "RECREATE");
              if (f.IsOpen())
              {
                if (hPhoRespSim)          hPhoRespSim->Write("h2_phoResp_truthVsReco");
                if (hPhoResp_measXtruth)  hPhoResp_measXtruth->Write("h2_phoResp_recoVsTruth");
                if (hPhoRecoData)         hPhoRecoData->Write("h_phoReco_data");
                if (hPhoRecoSim)          hPhoRecoSim->Write("h_phoReco_sim");
                if (hPhoTruthSim)         hPhoTruthSim->Write("h_phoTruth_sim");
                if (hPhoUnfoldTruth)      hPhoUnfoldTruth->Write("h_phoTruth_unfolded_data");
                if (hPhoUnfoldTruth_cov)  hPhoUnfoldTruth_cov->Write("h_phoTruth_unfolded_data_covariance");
                f.Close();
              }
            }
      
            // ----------------------------------------------------------------------
            // NEW (Step A validation): 1D photon unfolding QA package
            //
            // Output folder:
            //   <phoDir>/validation/
            //
            // Includes:
            //   (1) Iteration stability: rel(stat) + rel(change it vs it-1)
            //   (2) Closure: unfold(reco SIM) -> truth SIM
            //   (3) Half-closure: train(A) unfold(B) using binomial split of response
            // ----------------------------------------------------------------------
            {
              const string phoValDir = JoinPath(phoDir, "validation");
              EnsureDir(phoValDir);
      
                // -----------------------------
                // (1) Iteration stability (DATA)
                // -----------------------------
                if (hPhoRecoData && hPhoTruthSim && hPhoRecoSim && hPhoResp_measXtruth)
                {
                  const int kMaxIt = kMaxBayesIterPhoScan;
      
                  const vector<double>& xIt       = phoIterScan.xIt;
                  const vector<double>& exIt      = phoIterScan.exIt;
                  const vector<double>& yRelStat  = phoIterScan.yRelStat;
                  const vector<double>& eyRelStat = phoIterScan.eyRelStat;
                  const vector<double>& yRelDev   = phoIterScan.yRelChange;
                  const vector<double>& eyRelDev  = phoIterScan.eyRelChange;
                  const vector<double>& yQuad     = phoIterScan.yQuad;
                  const vector<double>& eyQuad    = phoIterScan.eyQuad;
                  const int bestIt                = phoIterScan.bestIt;
                  const double bestQuad           = phoIterScan.bestQuad;
      
                  if ((int)xIt.size() > 0)
                  {
                    TCanvas cSt("c_pho_iterStability_relChange_relStat", "c_pho_iterStability_relChange_relStat", 900, 700);
                    ApplyCanvasMargins1D(cSt);
      
                    // Determine y-range from content
                    double ymax = 0.0;
                    for (size_t i = 0; i < yRelStat.size(); ++i) ymax = std::max(ymax, yRelStat[i]);
                    for (size_t i = 0; i < yRelDev.size();  ++i) ymax = std::max(ymax, yRelDev[i]);
                    if (ymax <= 0.0) ymax = 0.1;
                    ymax *= 1.25;
      
                    TH1F frame("frame_phoIt", "", 1, 1.0, (double)kMaxIt + 0.5);
                    frame.SetMinimum(0.0);
                    frame.SetMaximum(ymax);
                    frame.SetTitle("");
                    frame.GetXaxis()->SetTitle("Iteration");
                    frame.GetYaxis()->SetTitle("Relative quantity");
                    frame.Draw("axis");
      
                    // total relative stat. uncertainty (RED)
                    TGraphErrors gStat((int)xIt.size(), &xIt[0], &yRelStat[0], &exIt[0], &eyRelStat[0]);
                    gStat.SetMarkerStyle(24);
                    gStat.SetMarkerSize(1.1);
                    gStat.SetMarkerColor(kRed + 1);
                    gStat.SetLineColor(kRed + 1);
                    gStat.SetLineWidth(2);
                    gStat.Draw("P same");
      
                    // total relative deviation (it vs it-1, with it=1 using 0->1 baseline) (BLUE)
                    TGraphErrors gDev((int)xIt.size(), &xIt[0], &yRelDev[0], &exIt[0], &eyRelDev[0]);
                    gDev.SetMarkerStyle(24);
                    gDev.SetMarkerSize(1.1);
                    gDev.SetMarkerColor(kBlue + 1);
                    gDev.SetLineColor(kBlue + 1);
                    gDev.SetLineWidth(2);
                    gDev.Draw("P same");
      
                    TLegend leg(0.55, 0.78, 0.89, 0.90);
                    leg.SetBorderSize(0);
                    leg.SetFillStyle(0);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.033);
                    leg.AddEntry(&gStat, "total relative stat. uncertainty", "p");
                    leg.AddEntry(&gDev,  "total relative deviation (it vs it-1)", "p");
                    leg.Draw();
      
                    {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(22);
                        tx.SetTextSize(0.040);
                        tx.DrawLatex(0.50, 0.965, "Iteration Stability, Photon 4 + MBD NS #geq 1, Run24pp");
                    }
      
                    // Label under legend (photon 1D unfolding)
                    {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(13);
                        tx.SetTextSize(0.032);
                        tx.DrawLatex(0.55, 0.74, "1D photon-yield unfolding");
                    }
      
                    SaveCanvas(cSt, JoinPath(phoValDir, "pho_unfold_iterStability_relChange_relStat.png"));
      
                    TCanvas cQuad("c_pho_iterStability_quadratureSum", "c_pho_iterStability_quadratureSum", 900, 700);
                    ApplyCanvasMargins1D(cQuad);
      
                    double ymaxQ = 0.0;
                    for (size_t i = 0; i < yQuad.size(); ++i) ymaxQ = std::max(ymaxQ, yQuad[i]);
                    if (ymaxQ <= 0.0) ymaxQ = 0.1;
                    ymaxQ *= 1.25;
      
                    TH1F frameQ("frame_phoItQuad", "", 1, 1.0, (double)kMaxIt + 0.5);
                    frameQ.SetMinimum(0.0);
                    frameQ.SetMaximum(ymaxQ);
                    frameQ.SetTitle("");
                    frameQ.GetXaxis()->SetTitle("Iteration");
                    frameQ.GetYaxis()->SetTitle("Quadrature sum");
                    frameQ.Draw("axis");
      
                    TGraphErrors gQuad((int)xIt.size(), &xIt[0], &yQuad[0], &exIt[0], &eyQuad[0]);
                    gQuad.SetMarkerStyle(20);
                    gQuad.SetMarkerSize(1.1);
                    gQuad.SetMarkerColor(kBlack);
                    gQuad.SetLineColor(kBlack);
                    gQuad.SetLineWidth(2);
                    gQuad.Draw("LP same");
      
                    TGraphErrors gBest;
                    if (bestIt > 0 && std::isfinite(bestQuad))
                    {
                      const double bestX[1]  = { (double)bestIt };
                      const double bestY[1]  = { bestQuad };
                      const double bestEX[1] = { 0.0 };
                      const double bestEY[1] = { 0.0 };
      
                      gBest = TGraphErrors(1, bestX, bestY, bestEX, bestEY);
                      gBest.SetMarkerStyle(29);
                      gBest.SetMarkerSize(1.6);
                      gBest.SetMarkerColor(kRed + 1);
                      gBest.SetLineColor(kRed + 1);
                      gBest.Draw("P same");
                    }
      
                    TLegend legQ(0.55, 0.78, 0.89, 0.90);
                    legQ.SetBorderSize(0);
                    legQ.SetFillStyle(0);
                    legQ.SetTextFont(42);
                    legQ.SetTextSize(0.033);
                    legQ.AddEntry(&gQuad, "quadrature sum", "lp");
                    if (bestIt > 0 && std::isfinite(bestQuad))
                    {
                      legQ.AddEntry(&gBest, TString::Format("minimum: it=%d", bestIt).Data(), "p");
                    }
                    legQ.Draw();
      
                    {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(22);
                        tx.SetTextSize(0.040);
                        tx.DrawLatex(0.50, 0.965, "Quadrature-sum optimization, Photon 4 + MBD NS #geq 1, Run24pp");
                    }
      
                    {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(13);
                        tx.SetTextSize(0.032);
                        tx.DrawLatex(0.55, 0.74, "1D photon-yield unfolding");
                    }
      
                    SaveCanvas(cQuad, JoinPath(phoValDir, "pho_unfold_iterStability_quadratureSum.png"));
      
                    vector<string> iterSummary;
                    iterSummary.push_back("Photon iteration-stability summary");
                    iterSummary.push_back(TString::Format("best iteration from quadrature sum = %d", bestIt).Data());
                    iterSummary.push_back(TString::Format("minimum quadrature sum = %.10g", bestQuad).Data());
                    iterSummary.push_back("relChange is defined as the successive-iterate difference it vs (it-1), with iteration 1 using the explicit 0->1 baseline.");
                    iterSummary.push_back("iteration 0 baseline is the measured input histogram mapped onto the unfolding comparison axis, so iteration 1 is included in the best-iteration choice.");
                    iterSummary.push_back("relStat is built from the full RooUnfold::kCovToy covariance restricted to the 9 analysis p_{T}^{#gamma} bins.");
                    iterSummary.push_back("support / underflow / overflow bins are excluded from the stability scalar.");
                    iterSummary.push_back("quadrature sum definition: sqrt(relStat^2 + relChange^2)");
                    iterSummary.push_back("");
      
                    for (size_t i = 0; i < xIt.size(); ++i)
                    {
                      iterSummary.push_back(
                        TString::Format("it=%d  relStat=%.10g  relChange=%.10g  quadratureSum=%.10g",
                                        (int)xIt[i], yRelStat[i], yRelDev[i], yQuad[i]).Data()
                      );
                    }
      
                    WriteTextFile(JoinPath(phoValDir, "pho_unfold_iterStability_bestIteration.txt"), iterSummary);
      
                    cout << ANSI_BOLD_CYN
                         << "[PHO ITER QA] best iteration from quadrature sum = " << bestIt
                         << "  minimum = " << bestQuad
                         << ANSI_RESET << "\n";
                  }
                }
      
                // -----------------------------
                // (2) Closure: unfold(reco SIM) -> truth SIM
                // -----------------------------
                if (hPhoRecoSim && hPhoTruthSim)
                {
                  RooUnfoldBayes uC(&respPho, hPhoRecoSim, kBayesIterPho);
                  uC.SetVerbose(0);
                  uC.SetNToys(kNToysPhoScan);
      
                  TH1* hUnfC = nullptr;
                  if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
                  hUnfC = uC.Hreco(RooUnfold::kCovToy);
                  if (gSystem) gSystem->RedirectOutput(0);
                  if (hUnfC)
                  {
                    hUnfC->SetDirectory(nullptr);
                    EnsureSumw2(hUnfC);
      
                    TH1* hRat = CloneTH1(hUnfC, "h_pho_closure_unfoldedOverTruth");
                    if (hRat)
                    {
                      hRat->SetDirectory(nullptr);
                      EnsureSumw2(hRat);
                      hRat->Divide(hPhoTruthSim);
      
                      hRat->SetTitle("");
                      hRat->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                      hRat->GetXaxis()->SetTitleOffset(1.15);
                      hRat->GetYaxis()->SetTitle("Closure: unfolded MC / truth MC");
                      hRat->SetMarkerStyle(24);
                      hRat->SetMarkerSize(1.1);
                      hRat->SetLineWidth(2);
        
                      TCanvas c("c_pho_closure", "c_pho_closure", 900, 700);
                      ApplyCanvasMargins1D(c);
      
                      const double xPlotMin = 10.0;
                      const double xPlotMax = 35.0;
      
                      const double ymin = 0.95;
                      const double ymax = 1.05;
      
                      hRat->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                      hRat->GetYaxis()->SetRangeUser(ymin, ymax);
                      hRat->Draw("E1");
      
                      TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");
      
                        // Top-left Bayes + Step label (match 2D closure styling)
                        {
                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextAlign(13);
                          tx.SetTextSize(0.038);
                          tx.DrawLatex(0.20, 0.90, TString::Format("Bayes it=%d (photon)", kBayesIterPho).Data());
                          tx.DrawLatex(0.20, 0.855, "1D photon-yield unfolding");
                         }
        
                          // Centered title
                          {
                            TLatex tx;
                            tx.SetNDC(true);
                            tx.SetTextFont(42);
                            tx.SetTextAlign(22);
                            tx.SetTextSize(0.040);
                            tx.DrawLatex(0.50, 0.965, "Photon closure test, unfold(reco) #rightarrow truth (SIM)");
                          }
        
                      SaveCanvas(c, JoinPath(phoValDir, "pho_closure_unfoldedOverTruth_vs_pTgamma.png"));
      
                      delete hRat;
                    }
      
                    delete hUnfC;
                  }
              }
      
                // -----------------------------
                // (3) Half-closure: train(A) unfold(B)
                //     Split response (reco x truth) bin-by-bin with Binomial(0.5)
                // -----------------------------
                if (hPhoResp_measXtruth && hPhoRecoSim && hPhoTruthSim)
                {
                  const int nReco  = hPhoResp_measXtruth->GetNbinsX();
                  const int nTruth = hPhoResp_measXtruth->GetNbinsY();
      
                  unsigned int seed = 1337u;
                  seed = seed * 131u + (unsigned int)nReco;
                  seed = seed * 131u + (unsigned int)nTruth;
                  std::mt19937 rng(seed);
      
                  TH2* hRspA = CloneTH2(hPhoResp_measXtruth, "h2_phoResp_recoVsTruth_halfA");
                  TH2* hRspB = CloneTH2(hPhoResp_measXtruth, "h2_phoResp_recoVsTruth_halfB");
      
                  TH1* hMeasA  = CloneTH1(hPhoRecoSim,  "h_phoRecoSim_halfA");
                  TH1* hMeasB  = CloneTH1(hPhoRecoSim,  "h_phoRecoSim_halfB");
                  TH1* hTruthA = CloneTH1(hPhoTruthSim, "h_phoTruthSim_halfA");
                  TH1* hTruthB = CloneTH1(hPhoTruthSim, "h_phoTruthSim_halfB");
      
                  if (hRspA)  { hRspA->SetDirectory(nullptr);  hRspA->Reset("ICES");  hRspA->Sumw2(); }
                  if (hRspB)  { hRspB->SetDirectory(nullptr);  hRspB->Reset("ICES");  hRspB->Sumw2(); }
                  if (hMeasA) { hMeasA->SetDirectory(nullptr); hMeasA->Reset("ICES"); EnsureSumw2(hMeasA); }
                  if (hMeasB) { hMeasB->SetDirectory(nullptr); hMeasB->Reset("ICES"); EnsureSumw2(hMeasB); }
                  if (hTruthA){ hTruthA->SetDirectory(nullptr);hTruthA->Reset("ICES");EnsureSumw2(hTruthA); }
                  if (hTruthB){ hTruthB->SetDirectory(nullptr);hTruthB->Reset("ICES");EnsureSumw2(hTruthB); }
      
                  if (!hRspA || !hRspB || !hMeasA || !hMeasB || !hTruthA || !hTruthB)
                  {
                    cout << ANSI_BOLD_YEL
                         << "[WARN] Photon half-closure: failed to allocate split response/marginals. Skipping.\n"
                         << ANSI_RESET;
                  }
                  else
                  {
                    vector<double> measA((std::size_t)nReco + 2, 0.0), measB((std::size_t)nReco + 2, 0.0);
                    vector<double> truA ((std::size_t)nTruth + 2, 0.0), truB ((std::size_t)nTruth + 2, 0.0);
      
                    for (int ix = 0; ix <= nReco + 1; ++ix)
                    {
                      for (int iy = 0; iy <= nTruth + 1; ++iy)
                      {
                        const double nRaw = hPhoResp_measXtruth->GetBinContent(ix, iy);
                        if (!(nRaw > 0.0))
                        {
                          hRspA->SetBinContent(ix, iy, 0.0); hRspA->SetBinError(ix, iy, 0.0);
                          hRspB->SetBinContent(ix, iy, 0.0); hRspB->SetBinError(ix, iy, 0.0);
                          continue;
                        }
      
                        const long long N = (long long)std::llround(nRaw);
                        if (N <= 0LL)
                        {
                          hRspA->SetBinContent(ix, iy, 0.0); hRspA->SetBinError(ix, iy, 0.0);
                          hRspB->SetBinContent(ix, iy, 0.0); hRspB->SetBinError(ix, iy, 0.0);
                          continue;
                        }
      
                        std::binomial_distribution<long long> d(N, 0.5);
                        const long long NA = d(rng);
                        const long long NB = N - NA;
      
                        hRspA->SetBinContent(ix, iy, (double)NA);
                        hRspA->SetBinError  (ix, iy, std::sqrt((double)NA));
      
                        hRspB->SetBinContent(ix, iy, (double)NB);
                        hRspB->SetBinError  (ix, iy, std::sqrt((double)NB));
      
                        if (ix >= 0 && ix <= nReco + 1)
                        {
                          measA[(std::size_t)ix] += (double)NA;
                          measB[(std::size_t)ix] += (double)NB;
                        }
                        if (iy >= 0 && iy <= nTruth + 1)
                        {
                          truA[(std::size_t)iy] += (double)NA;
                          truB[(std::size_t)iy] += (double)NB;
                        }
                      }
                    }
      
                    for (int ib = 0; ib <= nReco + 1; ++ib)
                    {
                      hMeasA->SetBinContent(ib, measA[(std::size_t)ib]);
                      hMeasA->SetBinError  (ib, std::sqrt(measA[(std::size_t)ib]));
                      hMeasB->SetBinContent(ib, measB[(std::size_t)ib]);
                      hMeasB->SetBinError  (ib, std::sqrt(measB[(std::size_t)ib]));
                    }
      
                    for (int ib = 0; ib <= nTruth + 1; ++ib)
                    {
                      hTruthA->SetBinContent(ib, truA[(std::size_t)ib]);
                      hTruthA->SetBinError  (ib, std::sqrt(truA[(std::size_t)ib]));
                      hTruthB->SetBinContent(ib, truB[(std::size_t)ib]);
                      hTruthB->SetBinError  (ib, std::sqrt(truB[(std::size_t)ib]));
                    }
      
                    RooUnfoldResponse respPhoA(hMeasA, hTruthA, hRspA, "respPho_halfA", "respPho_halfA");
      
                    RooUnfoldBayes uH(&respPhoA, hMeasB, kBayesIterPho);
                    uH.SetVerbose(0);
                    uH.SetNToys(kNToysPhoScan);
      
                    TH1* hUnfB = nullptr;
                    if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
                    hUnfB = uH.Hreco(RooUnfold::kCovToy);
                    if (gSystem) gSystem->RedirectOutput(0);
                    if (hUnfB)
                    {
                      hUnfB->SetDirectory(nullptr);
                      EnsureSumw2(hUnfB);
      
                      TH1* hRat = CloneTH1(hUnfB, "h_pho_halfClosure_unfoldedOverTruth");
                      if (hRat)
                      {
                        hRat->SetDirectory(nullptr);
                        EnsureSumw2(hRat);
                        hRat->Divide(hTruthB);
      
                        hRat->SetTitle("");
                        hRat->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hRat->GetXaxis()->SetTitleOffset(1.15);
                        hRat->GetYaxis()->SetTitle("Half-closure: unfolded MC / truth MC");
                        hRat->SetMarkerStyle(24);
                        hRat->SetMarkerSize(1.1);
                        hRat->SetLineWidth(2);
        
                        TCanvas c("c_pho_halfClosure", "c_pho_halfClosure", 900, 700);
                        ApplyCanvasMargins1D(c);
      
                        const double xPlotMin = 10.0;
                        const double xPlotMax = 35.0;
      
                        double ymin =  1e99;
                        double ymax = -1e99;
                        bool havePoint = false;
      
                        for (int ib = 1; ib <= hRat->GetNbinsX(); ++ib)
                        {
                          const double lo = hRat->GetXaxis()->GetBinLowEdge(ib);
                          const double hi = hRat->GetXaxis()->GetBinUpEdge(ib);
                          if (lo < xPlotMin || hi > xPlotMax) continue;
      
                          const double y  = hRat->GetBinContent(ib);
                          const double ey = hRat->GetBinError(ib);
                          if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                          if (y <= 0.0 && ey <= 0.0) continue;
      
                          havePoint = true;
                          ymin = std::min(ymin, y - ey);
                          ymax = std::max(ymax, y + ey);
                        }
      
                        if (!havePoint || !(ymin < ymax))
                        {
                          ymin = 0.90;
                          ymax = 1.10;
                        }
                        else
                        {
                          const double pad = std::max(0.15 * (ymax - ymin), 0.01);
                          ymin -= pad;
                          ymax += pad;
      
                          if (ymin > 1.0) ymin = 1.0 - pad;
                          if (ymax < 1.0) ymax = 1.0 + pad;
                        }
      
                        hRat->GetXaxis()->SetRangeUser(xPlotMin, xPlotMax);
                        hRat->GetYaxis()->SetRangeUser(ymin, ymax);
                        hRat->Draw("E1");
      
                        TLine l1(xPlotMin, 1.0, xPlotMax, 1.0);
                        l1.SetLineStyle(2);
                        l1.SetLineWidth(2);
                        l1.Draw("same");
      
                        // Top-left Bayes + Step label (match 2D half-closure styling)
                        {
                            TLatex tx;
                            tx.SetNDC();
                            tx.SetTextFont(42);
                            tx.SetTextAlign(13);
                            tx.SetTextSize(0.038);
                            tx.DrawLatex(0.20, 0.90, TString::Format("Bayes it=%d (photon)", kBayesIterPho).Data());
                            tx.DrawLatex(0.20, 0.855, "1D photon-yield unfolding");
                        }
      
                        // Centered title
                        {
                            TLatex tx;
                            tx.SetNDC(true);
                            tx.SetTextFont(42);
                            tx.SetTextAlign(22);
                            tx.SetTextSize(0.040);
                            tx.DrawLatex(0.50, 0.965, "Photon half-closure test, train(A) unfold(B) (SIM)");
                        }
      
                        SaveCanvas(c, JoinPath(phoValDir, "pho_halfClosure_unfoldedOverTruth_vs_pTgamma.png"));
      
                        delete hRat;
                      }
      
                      delete hUnfB;
                    }
                  }
      
                  if (hTruthB) delete hTruthB;
                  if (hTruthA) delete hTruthA;
                  if (hMeasB)  delete hMeasB;
                  if (hMeasA)  delete hMeasA;
                  if (hRspB)   delete hRspB;
                  if (hRspA)   delete hRspA;
                }
            }

            {
              const string srcLog = JoinPath(phoDir, "pho_unfolded_truth_pTgamma_overlay_logy.png");
              const string dstLog = JoinPath(phoYieldDir, "pho_unfolded_truth_pTgamma_overlay_logy.png");

              if (gSystem)
              {
                gSystem->Unlink(dstLog.c_str());
                gSystem->Rename(srcLog.c_str(), dstLog.c_str());
              }
            }

            delete hPhoRecoData;
            delete hPhoRecoSim;
            delete hPhoTruthSim;
            delete hPhoRespSim;
            delete hPhoResp_measXtruth;
            if (hPhoUnfoldTruth) delete hPhoUnfoldTruth;
            if (hPhoUnfoldTruth_cov) delete hPhoUnfoldTruth_cov;
          }
          while (false);
        }

        // -------------------------------------------------------------------------
        // (B) Jet unfolding per radius: unfold global-bin vector and unflatten back
        // -------------------------------------------------------------------------
        const int kDefaultBayesIterXJ = 3;
        const int kMaxBayesIterXJScan = 10;

        // Toy settings for xJ:
        //   - "final" is for the baseline unfolded spectrum that feeds your main outputs
        //   - "scan" is for iteration-stability / closure / half-closure diagnostics
        const int kNToysXJFinal = 600;
        const int kNToysXJScan  = 120;

      std::map<std::string, std::vector<TH1*>> perPhoHistsByRKey;

        for (const auto& rKey : kRKeys)
        {
          const double R = RFromKey(rKey);
          const string rOut = JoinPath(outBase, rKey);
          EnsureDir(rOut);

            const string overlayOut = JoinPath(rOut, "LHC_overlay");
            EnsureDir(overlayOut);

            const string toyVsCovOut = JoinPath(rOut, "ToyUnfoldingVsCovariance");
            EnsureDir(toyVsCovOut);

            const string beforeAfterDataOut  = JoinPath(rOut, "before_after_unfoldingOverlay_data");
            const string beforeAfterTruthOut = JoinPath(rOut, "before_after_unfoldingOverlay_truth");
            const string withAndWithoutJetEffOut = JoinPath(rOut, "withAndWithoutJetEffCorr");
            const string jetEffQAOut = JoinPath(withAndWithoutJetEffOut, "jetEffQA");
            const string finalUnfoldedOut = JoinPath(rOut, "finalUnfolded_dNdXJ");
            const string closurePlotsOut = JoinPath(rOut, "closurePlots");
            const string iterationStabilityOut = JoinPath(rOut, "iterationStabilityCheck");
            const string purityCorrQAOut = JoinPath(rOut, "purityCorrQA");
            EnsureDir(beforeAfterDataOut);
            EnsureDir(beforeAfterTruthOut);
            EnsureDir(finalUnfoldedOut);
            EnsureDir(closurePlotsOut);
            EnsureDir(iterationStabilityOut);
            if (gApplyPurityCorrectionForUnfolding)
            {
              EnsureDir(withAndWithoutJetEffOut);
              EnsureDir(jetEffQAOut);
              EnsureDir(purityCorrQAOut);
            }

            auto MakePuritySubtractionQA =
              [&](TH2* h2SideC)->void
            {
              if (!gApplyPurityCorrectionForUnfolding || !h2SideC) return;

              const auto& analysisRecoBins = UnfoldAnalysisRecoPtBins();
              const int nPtPads = (int)analysisRecoBins.size();
              const int nPtCols = 3;
              const int nPtRows = 3;

              auto DrawPuritySubtractionPad =
                [&](int iPad)->TH1D*
              {
                if (iPad < 0 || iPad >= nPtPads) return nullptr;

                const PtBin& b = analysisRecoBins[iPad];
                const double cen = 0.5 * (b.lo + b.hi);
                const int ix = h2SideC->GetXaxis()->FindBin(cen);

                TH1D* hPerSideC = nullptr;
                if (ix >= 1 && ix <= h2SideC->GetNbinsX())
                {
                  hPerSideC = h2SideC->ProjectionY(
                    TString::Format("h_puritySubtractionInput_%s_pTbin%d", rKey.c_str(), iPad + 1).Data(),
                    ix, ix, "e"
                  );
                  if (hPerSideC)
                  {
                    hPerSideC->SetDirectory(nullptr);
                    EnsureSumw2(hPerSideC);
                  }
                }

                if (hPerSideC)
                {
                  double maxY = 0.0;
                  const int nxb = hPerSideC->GetNbinsX();
                  for (int ib = 1; ib <= nxb; ++ib)
                  {
                    const double y  = hPerSideC->GetBinContent(ib);
                    const double ey = hPerSideC->GetBinError(ib);
                    const double v  = y + ey;
                    if (v > maxY) maxY = v;
                  }

                  hPerSideC->SetMinimum(0.0);
                  hPerSideC->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                  hPerSideC->GetXaxis()->SetRangeUser(0.0, 2.0);
                  hPerSideC->SetTitle("");
                  hPerSideC->GetXaxis()->SetTitle("x_{J}");
                  hPerSideC->GetYaxis()->SetTitle("Counts");
                  hPerSideC->Draw("E1");
                }
                else
                {
                  TH1F frame("frame_puritySubtractionInput","", 1, 0.0, 2.0);
                  frame.SetMinimum(0.0);
                  frame.SetMaximum(1.0);
                  frame.SetTitle("");
                  frame.GetXaxis()->SetTitle("x_{J}");
                  frame.GetYaxis()->SetTitle("Counts");
                  frame.Draw("axis");

                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextSize(0.050);
                  tx.DrawLatex(0.16, 0.50, "MISSING");
                }

                {
                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(22);
                  tx.SetTextSize(0.042);

                  tx.DrawLatex(0.52, 0.955,
                               TString::Format("Purity-subtraction input x_{J#gamma} (region C), p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                               b.lo, b.hi, R).Data());
                }

                {
                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(13);
                  tx.SetTextSize(0.040);

                  const double xL = 0.14;
                  tx.DrawLatex(xL, 0.88, "Run24pp");
                  tx.DrawLatex(xL, 0.81, "Trigger = Photon 4 + MBD NS #geq 1");
                  tx.DrawLatex(xL, 0.74, "#Delta #phi > 7#pi/8");
                  tx.DrawLatex(xL, 0.67, "z_{vtx} < 60 cm");
                  tx.DrawLatex(xL, 0.60, "Region C = iso + non-tight");
                }

                return hPerSideC;
              };

              TCanvas cTbl(
                TString::Format("c_puritySubtractionInput_tbl_%s", rKey.c_str()).Data(),
                "c_puritySubtractionInput_tbl",
                1800, 1300
              );
              cTbl.Divide(nPtCols, nPtRows, 0.001, 0.001);

              std::vector<TH1*> keepPurityTable;
              keepPurityTable.reserve((std::size_t)nPtPads);

              for (int ipad = 0; ipad < nPtPads; ++ipad)
              {
                cTbl.cd(ipad + 1);
                gPad->SetLeftMargin(0.12);
                gPad->SetRightMargin(0.04);
                gPad->SetBottomMargin(0.12);
                gPad->SetTopMargin(0.06);

                TH1D* hKeep = DrawPuritySubtractionPad(ipad);
                if (hKeep) keepPurityTable.push_back(hKeep);
              }

              SaveCanvas(cTbl, JoinPath(purityCorrQAOut, "table3x3_puritySubtractionInput_xJ.png"));

              for (auto* h : keepPurityTable) if (h) delete h;
              keepPurityTable.clear();

              for (int i = 0; i < nPtPads; ++i)
              {
                TCanvas cSingle(
                  TString::Format("c_puritySubtractionInput_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
                  "c_puritySubtractionInput_single",
                  900, 700
                );
                ApplyCanvasMargins1D(cSingle);

                TH1D* hKeep = DrawPuritySubtractionPad(i);

                SaveCanvas(
                  cSingle,
                  JoinPath(
                    purityCorrQAOut,
                    TString::Format("puritySubtractionInput_xJ_pTbin%d.png", i + 1).Data()
                  )
                );

                if (hKeep) delete hKeep;
              }
            };

            const string nameReco      = "h2_unfoldReco_pTgamma_xJ_incl_"           + rKey;
            const string nameRecoC     = "h2_unfoldReco_pTgamma_xJ_incl_sidebandC_" + rKey;
            const string nameTruth     = "h2_unfoldTruth_pTgamma_xJ_incl_"          + rKey;
            const string nameRsp       = "h2_unfoldResponse_pTgamma_xJ_incl_"       + rKey;
            const string nameJetEffDen = "h2_unfoldJetEffDen_pTgamma_xJ_incl_"      + rKey;
            const string nameJetEffNum = "h2_unfoldJetEffNum_pTgamma_xJ_incl_"      + rKey;

            TH2* h2RecoData_in       = GetObj<TH2>(dsData, nameReco,  true, true, true);
            TH2* h2RecoData_sideC_in = (gApplyPurityCorrectionForUnfolding
                                        ? GetObj<TH2>(dsData, nameRecoC, true, true, true)
                                        : nullptr);
            TH2* h2RecoSim_in        = GetObj<TH2>(dsSim,  nameReco,      true,  true, true);
            TH2* h2TruthSim_in       = GetObj<TH2>(dsSim,  nameTruth,     true,  true, true);
            TH2* h2RspSim_in         = GetObj<TH2>(dsSim,  nameRsp,       true,  true, true);
            TH2* h2JetEffDen_in      = (gApplyPurityCorrectionForUnfolding
                                        ? GetObj<TH2>(dsSim, nameJetEffDen, true, true, false)
                                        : nullptr);
            TH2* h2JetEffNum_in      = (gApplyPurityCorrectionForUnfolding
                                        ? GetObj<TH2>(dsSim, nameJetEffNum, true, true, false)
                                        : nullptr);

          if (!h2RecoData_in || !h2RecoSim_in || !h2TruthSim_in || !h2RspSim_in ||
              (gApplyPurityCorrectionForUnfolding && !h2RecoData_sideC_in))
          {
              auto Status = [&](Dataset& ds, const string& relName, TObject* obj)->string
              {
                if (obj) return "FOUND";
                const string fp = FullPath(ds, relName);
                auto it = ds.missingReason.find(fp);
                if (it == ds.missingReason.end()) return "MISSING";
                return string("MISSING (") + it->second + ")";
              };

              cout << ANSI_BOLD_YEL
                   << "[WARN] Skipping RooUnfold xJ for " << rKey << " due to missing/zero inputs:\n"
                   << "       DATA  " << FullPath(dsData, nameReco)  << " : " << Status(dsData, nameReco,  h2RecoData_in)  << "\n";
              if (gApplyPurityCorrectionForUnfolding)
                cout << "       DATA  " << FullPath(dsData, nameRecoC) << " : " << Status(dsData, nameRecoC, h2RecoData_sideC_in) << "\n";
              cout << "       SIM   " << FullPath(dsSim,  nameReco)  << " : " << Status(dsSim,  nameReco,  h2RecoSim_in)   << "\n"
                   << "       SIM   " << FullPath(dsSim,  nameTruth) << " : " << Status(dsSim,  nameTruth, h2TruthSim_in)  << "\n"
                   << "       SIM   " << FullPath(dsSim,  nameRsp)   << " : " << Status(dsSim,  nameRsp,   h2RspSim_in)    << "\n"
                   << ANSI_RESET << "\n";
              continue;
          }

          TH2* h2RecoData  = CloneTH2(h2RecoData_in,  TString::Format("h2RecoData_%s", rKey.c_str()).Data());
          TH2* h2RecoSim   = CloneTH2(h2RecoSim_in,   TString::Format("h2RecoSim_%s",  rKey.c_str()).Data());
          TH2* h2TruthSim  = CloneTH2(h2TruthSim_in,  TString::Format("h2TruthSim_%s", rKey.c_str()).Data());
          TH2* h2RspSim    = CloneTH2(h2RspSim_in,    TString::Format("h2RspSim_truthVsReco_%s", rKey.c_str()).Data());

          EnsureSumw2(h2RecoData);
          EnsureSumw2(h2RecoSim);
          EnsureSumw2(h2TruthSim);
          EnsureSumw2(h2RspSim);

          if (gApplyPurityCorrectionForUnfolding)
          {
            cout << ANSI_BOLD_CYN
                 << "\n[DEBUG PURITY XJ] Building purity-corrected reco input for " << rKey << "\n"
                 << "  source reco hist   = " << (h2RecoData ? h2RecoData->GetName() : "<null>") << "\n"
                 << "  source sideband-C  = " << (h2RecoData_sideC_in ? h2RecoData_sideC_in->GetName() : "<null>") << "\n"
                 << ANSI_RESET;

            DumpTH2Summary(TString::Format("DATA reco BEFORE purity correction (%s)", rKey.c_str()).Data(), h2RecoData);
            DumpTH2Summary(TString::Format("DATA sideband-C input (%s)", rKey.c_str()).Data(), h2RecoData_sideC_in);
            MakePuritySubtractionQA(h2RecoData_sideC_in);

            for (int iDbg = 0; iDbg < kNPtBins; ++iDbg)
            {
              double ADbg = 0.0, BDbg = 0.0, CDbg = 0.0, DDbg = 0.0;
              double SADbg = 0.0, eSADbg = 0.0, nbkgADbg = 0.0;
              ComputeABCDSignalCounts(iDbg, ADbg, BDbg, CDbg, DDbg, SADbg, eSADbg, nbkgADbg);

              const PtBin& bDbg = PtBins()[iDbg];
              const double cenDbg = 0.5 * (bDbg.lo + bDbg.hi);
              const int ixDbg = h2RecoData ? h2RecoData->GetXaxis()->FindBin(cenDbg) : -1;

              cout << "  [ABCD " << rKey << "] pT=" << bDbg.lo << "-" << bDbg.hi
                   << "  cen=" << cenDbg
                   << "  ixReco=" << ixDbg
                   << "  A=" << ADbg
                   << "  B=" << BDbg
                   << "  C=" << CDbg
                   << "  D=" << DDbg
                   << "  SA=" << SADbg
                   << "  eSA=" << eSADbg
                   << "  nbkgA=" << nbkgADbg
                   << "  scaleC=" << ((CDbg > 0.0) ? (nbkgADbg / CDbg) : 0.0)
                   << "\n";
            }

            if (!ApplyPurityCorrectionToRecoXJHist(h2RecoData, h2RecoData_sideC_in))
            {
              cout << ANSI_BOLD_RED
                   << "[ERROR] Failed to build purity-corrected xJ reco input for " << rKey << ". Skipping this radius."
                   << ANSI_RESET << "\n";
              delete h2RecoData;
              delete h2RecoSim;
              delete h2TruthSim;
              delete h2RspSim;
              continue;
            }

            DumpTH2Summary(TString::Format("DATA reco AFTER purity correction (%s)", rKey.c_str()).Data(), h2RecoData);

            for (int iDbg = 0; iDbg < kNPtBins; ++iDbg)
            {
              const PtBin& bDbg = PtBins()[iDbg];
              const double cenDbg = 0.5 * (bDbg.lo + bDbg.hi);
              const int ixDbg = h2RecoData ? h2RecoData->GetXaxis()->FindBin(cenDbg) : -1;
              DumpProjectionYSummary(
                TString::Format("Purity-corrected reco xJ projection (%s, pT %d-%d)", rKey.c_str(), bDbg.lo, bDbg.hi).Data(),
                h2RecoData,
                ixDbg
              );
            }
          }
          else
          {
              cout << ANSI_BOLD_CYN
                   << "\n[DEBUG XJ INPUT] Using NON-purity-corrected reco input for " << rKey << "\n"
                   << ANSI_RESET;
              DumpTH2Summary(TString::Format("DATA reco input without purity correction (%s)", rKey.c_str()).Data(), h2RecoData);
            }

            // --- Optional: ATLAS-style combinatoric jet background subtraction ---
            // Subtract H_comb^{embed} from the reco data histogram before unfolding.
            // The combinatoric template lives in the SIM (embedded) file.
            if (gApplyCombinatoricSubtractionForUnfolding && IsEmbeddedSimSample(CurrentSimSample()))
            {
              const string nameComb = "h2_unfoldRecoCombinatoric_pTgamma_xJ_incl_" + rKey;
              TH2* h2Comb_in = GetObj<TH2>(dsSim, nameComb, true, true, false);

              if (h2Comb_in)
              {
                TH2* h2Comb = CloneTH2(h2Comb_in, TString::Format("h2Comb_%s", rKey.c_str()).Data());
                EnsureSumw2(h2Comb);

                DumpTH2Summary(TString::Format("SIM combinatoric template BEFORE subtraction (%s)", rKey.c_str()).Data(), h2Comb);
                DumpTH2Summary(TString::Format("DATA reco BEFORE combinatoric subtraction (%s)", rKey.c_str()).Data(), h2RecoData);

                h2RecoData->Add(h2Comb, -1.0);

                DumpTH2Summary(TString::Format("DATA reco AFTER combinatoric subtraction (%s)", rKey.c_str()).Data(), h2RecoData);

                cout << ANSI_BOLD_CYN
                     << "[COMB SUB] Subtracted combinatoric template from reco data for " << rKey << "\n"
                     << ANSI_RESET;

                delete h2Comb;
              }
              else
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] gApplyCombinatoricSubtractionForUnfolding=true but combinatoric histogram not found: "
                     << nameComb << " in SIM file. Proceeding without combinatoric subtraction for " << rKey << ".\n"
                     << ANSI_RESET;
              }
            }

        DumpTH2Summary(TString::Format("SIM reco input (%s)", rKey.c_str()).Data(), h2RecoSim);
        DumpTH2Summary(TString::Format("SIM truth input (%s)", rKey.c_str()).Data(), h2TruthSim);
        DumpTH2Summary(TString::Format("SIM response truthVsReco (%s)", rKey.c_str()).Data(), h2RspSim);

        const int nGlobTruth_rsp = h2RspSim->GetNbinsX();
        const int nGlobReco_rsp  = h2RspSim->GetNbinsY();

        TH1D* hMeasSimGlob  = FlattenTH2ToGlobal(h2RecoSim,  TString::Format("hMeasSimGlob_%s",  rKey.c_str()).Data());
        TH1D* hTruthSimGlob = FlattenTH2ToGlobal(h2TruthSim, TString::Format("hTruthSimGlob_%s", rKey.c_str()).Data());
        TH1D* hMeasDataGlob = FlattenTH2ToGlobal(h2RecoData, TString::Format("hMeasDataGlob_%s", rKey.c_str()).Data());

        if (!hMeasSimGlob || !hTruthSimGlob || !hMeasDataGlob)
        {
          cout << ANSI_BOLD_YEL << "[WARN] Failed to flatten 2D hists for " << rKey << ". Skipping." << ANSI_RESET << "\n";
          if (hMeasSimGlob)  delete hMeasSimGlob;
          if (hTruthSimGlob) delete hTruthSimGlob;
          if (hMeasDataGlob) delete hMeasDataGlob;
          delete h2RecoData;
          delete h2RecoSim;
          delete h2TruthSim;
          delete h2RspSim;
          continue;
        }

        DumpTH1Summary(TString::Format("Flattened SIM reco global (%s)", rKey.c_str()).Data(), hMeasSimGlob);
        DumpTH1Summary(TString::Format("Flattened SIM truth global (%s)", rKey.c_str()).Data(), hTruthSimGlob);
        DumpTH1Summary(TString::Format("Flattened DATA reco global (%s)", rKey.c_str()).Data(), hMeasDataGlob);

        if (hMeasSimGlob->GetNbinsX() != nGlobReco_rsp || hTruthSimGlob->GetNbinsX() != nGlobTruth_rsp)
        {
          cout << ANSI_BOLD_RED
               << "[ERROR] Global-bin size mismatch for " << rKey << ":\n"
               << "        response nGlobTruth=" << nGlobTruth_rsp << " nGlobReco=" << nGlobReco_rsp << "\n"
               << "        flattened SIM truth nbins=" << hTruthSimGlob->GetNbinsX() << "\n"
               << "        flattened SIM reco  nbins=" << hMeasSimGlob ->GetNbinsX() << "\n"
               << "        Aborting this radius."
               << ANSI_RESET << "\n";

          delete hMeasSimGlob;
          delete hTruthSimGlob;
          delete hMeasDataGlob;
          delete h2RecoData;
          delete h2RecoSim;
          delete h2TruthSim;
          delete h2RspSim;
          continue;
        }

        TH2D* hRsp_measXtruth = new TH2D(
          TString::Format("h2_unfoldResponse_global_recoVsTruth_%s", rKey.c_str()).Data(),
          "Response matrix; global bin (reco); global bin (truth)",
          nGlobReco_rsp,  -0.5, nGlobReco_rsp  - 0.5,
          nGlobTruth_rsp, -0.5, nGlobTruth_rsp - 0.5
        );
        hRsp_measXtruth->SetDirectory(nullptr);
        hRsp_measXtruth->Sumw2();

        for (int ixTruth = 0; ixTruth <= nGlobTruth_rsp + 1; ++ixTruth)
        {
          for (int iyReco = 0; iyReco <= nGlobReco_rsp + 1; ++iyReco)
          {
            hRsp_measXtruth->SetBinContent(iyReco, ixTruth, h2RspSim->GetBinContent(ixTruth, iyReco));
            hRsp_measXtruth->SetBinError  (iyReco, ixTruth, h2RspSim->GetBinError  (ixTruth, iyReco));
          }
        }

          RooUnfoldResponse respXJ(hMeasSimGlob, hTruthSimGlob, hRsp_measXtruth,
                                  TString::Format("respXJ_%s", rKey.c_str()).Data(),
                                  TString::Format("respXJ_%s", rKey.c_str()).Data());

          const IterScanSummary xjIterScan = BuildXJIterScan(
            respXJ,
            hMeasDataGlob,
            hTruthSimGlob,
            h2RecoData,
            h2TruthSim,
            kMaxBayesIterXJScan,
            kNToysXJScan,
            rKey
          );

          const int kBayesIterXJ = (xjIterScan.bestIt > 0 ? xjIterScan.bestIt : kDefaultBayesIterXJ);

          cout << ANSI_BOLD_CYN
               << "[UNF ITER AUTO] rKey=" << rKey
               << "  Using Bayes iteration " << kBayesIterXJ;
          if (xjIterScan.bestIt > 0 && std::isfinite(xjIterScan.bestQuad))
          {
            cout << " from quadrature-sum minimum (" << xjIterScan.bestQuad << ")";
          }
          else
          {
            cout << " (fallback default; automatic scan unavailable)";
          }
          cout << ANSI_RESET << "\n";

          RooUnfoldBayes unfoldXJ_toy(&respXJ, hMeasDataGlob, kBayesIterXJ);
          unfoldXJ_toy.SetVerbose(0);
          unfoldXJ_toy.SetNToys(kNToysXJFinal);

          TH1* hUnfoldTruthGlob = nullptr;
          if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
          hUnfoldTruthGlob = unfoldXJ_toy.Hreco(RooUnfold::kCovToy);
          if (gSystem) gSystem->RedirectOutput(0);
          if (hUnfoldTruthGlob) hUnfoldTruthGlob->SetDirectory(nullptr);

          RooUnfoldBayes unfoldXJ_cov(&respXJ, hMeasDataGlob, kBayesIterXJ);
          unfoldXJ_cov.SetVerbose(0);

          TH1* hUnfoldTruthGlob_cov = unfoldXJ_cov.Hreco(RooUnfold::kCovariance);
          if (hUnfoldTruthGlob_cov) hUnfoldTruthGlob_cov->SetDirectory(nullptr);

        if (!hUnfoldTruthGlob)
          {
            cout << ANSI_BOLD_RED << "[ERROR] RooUnfold returned nullptr truth histogram for " << rKey << ". Skipping." << ANSI_RESET << "\n";
            delete hMeasSimGlob;
            delete hTruthSimGlob;
            delete hMeasDataGlob;
            delete hRsp_measXtruth;
            delete h2RecoData;
            delete h2RecoSim;
            delete h2TruthSim;
            delete h2RspSim;
            continue;
          }

          TH2* h2UnfoldTruth = UnflattenGlobalToTH2(
            hUnfoldTruthGlob,
            h2TruthSim,
            TString::Format("h2_unfoldedTruth_pTgamma_xJ_incl_%s", rKey.c_str()).Data()
          );

          TH2* h2UnfoldTruth_cov = nullptr;
          if (hUnfoldTruthGlob_cov)
          {
            h2UnfoldTruth_cov = UnflattenGlobalToTH2(
              hUnfoldTruthGlob_cov,
              h2TruthSim,
              TString::Format("h2_unfoldedTruth_pTgamma_xJ_incl_cov_%s", rKey.c_str()).Data()
            );
        }

        // Save response scatter (SIM) into the DATA unfolding folder for convenience
        {
          TCanvas c(TString::Format("c_rsp_scatter_%s", rKey.c_str()).Data(), "c_rsp_scatter", 900, 800);
          ApplyCanvasMargins2D(c);

          TH2* hScat = CloneTH2(h2RspSim, TString::Format("h2RspScat_%s", rKey.c_str()).Data());
          if (hScat)
          {
            c.SetLogz(1);

              hScat->SetTitle("");
              hScat->GetXaxis()->SetTitle("global bin (truth: p_{T}^{#gamma}, x_{J})");
              hScat->GetYaxis()->SetTitle("global bin (reco:  p_{T}^{#gamma}, x_{J})");
              hScat->GetZaxis()->SetTitle("Counts");

              hScat->GetXaxis()->SetTitleOffset(1.35);
              hScat->GetYaxis()->SetTitleOffset(1.45);

              hScat->GetZaxis()->SetTitleSize(0.035);
              hScat->GetZaxis()->SetTitleOffset(1.25);
              hScat->GetZaxis()->SetLabelSize(0.030);

              // required for log-z (must be > 0)
              hScat->SetMinimum(0.5);

              hScat->Draw("COLZ");

              {
                  TString simMergeLabel;
                  if      (allPhoton5and10and20sim) simMergeLabel = "Photon 5+10+20 GeV";
                  else if (bothPhoton5and10sim)     simMergeLabel = "Photon 5+10 GeV";
                  else if (bothPhoton5and20sim)     simMergeLabel = "Photon 5+20 GeV";
                  else if (bothPhoton10and20sim)    simMergeLabel = "Photon 10+20 GeV";
                  else                              simMergeLabel = "Photon";

                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(22);
                  tx.SetTextSize(0.040);
                  tx.DrawLatex(0.50, 0.965,
                               TString::Format("2D Response Matrix, %s, R = %.1f", simMergeLabel.Data(), R).Data());
              }

            SaveCanvas(c, JoinPath(rOut, "unfold_response_globalTruth_vs_globalReco_SCAT.png"));
            delete hScat;
          }
        }

        // Save unfolded 2D truth map
        if (h2UnfoldTruth)
        {
          TCanvas c(TString::Format("c_unf2D_%s", rKey.c_str()).Data(), "c_unf2D", 950, 800);
          ApplyCanvasMargins2D(c);
          c.SetLogz();

          h2UnfoldTruth->SetTitle("");
          h2UnfoldTruth->GetXaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
          h2UnfoldTruth->GetYaxis()->SetTitle("x_{J}^{truth}");
          h2UnfoldTruth->Draw("colz");

          DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsData), 0.034, 0.045);
          DrawLatexLines(0.14,0.78, { TString::Format("Unfolded truth (particle-level) (%s, R=%.1f)", rKey.c_str(), R).Data(),
                                      TString::Format("Bayes it=%d", kBayesIterXJ).Data() }, 0.030, 0.040);

          SaveCanvas(c, JoinPath(rOut, "unfoldedTruth_pTgamma_xJ_colz.png"));
        }

        vector<string> lines;
        lines.push_back("RooUnfold pipeline summary");
        lines.push_back(TString::Format("Radius: %s (R=%.1f)", rKey.c_str(), R).Data());
        lines.push_back(TString::Format("Photon unfolding: Bayes it=%d (kCovToy, Ntoys=%d) [also kCovariance for comparisons]", kBayesIterPho, kNToysPhoFinal).Data());
        lines.push_back(TString::Format("xJ unfolding (global-bin): Bayes it=%d (kCovToy, Ntoys=%d) [also kCovariance for comparisons]", kBayesIterXJ, kNToysXJFinal).Data());
        lines.push_back("");

        const auto& analysisRecoBins = UnfoldAnalysisRecoPtBins();
        const int nPtAll = (int)analysisRecoBins.size();
        const int nPtCols = 3;
        const int nPtRows = 3;
        const int nPtPads = nPtCols * nPtRows;

          vector<TH1*> perPhoHists(nPtAll, nullptr);
          vector<TH1*> perPhoHists_cov(nPtAll, nullptr);
          vector<TH1*> perPhoErrRatio(nPtAll, nullptr);
          vector<TH1*> perPhoBeforeDataHists(nPtAll, nullptr);
          vector<TH1*> perPhoTruthHists(nPtAll, nullptr);
          vector<TH1*> perPhoHists_jetEffCorr(nPtAll, nullptr);
          vector<TH1*> ratioBeforeVsAfterHists(nPtAll, nullptr);
          vector<TH1*> ratioTruthVsUnfoldedHists(nPtAll, nullptr);

          TH2* h2JetEff = nullptr;
          TH2* h2UnfoldTruth_jetEffCorr = nullptr;

          if (gApplyPurityCorrectionForUnfolding)
          {
            if (!h2JetEffDen_in || !h2JetEffNum_in)
            {
              cout << ANSI_BOLD_YEL
                   << "[WARN] Missing jet-efficiency inputs for " << rKey
                   << " (need " << nameJetEffDen << " and " << nameJetEffNum << "). "
                   << "Skipping withAndWithoutJetEffCorr outputs for this radius."
                   << ANSI_RESET << "\n";
            }
            else
            {
                auto DrawJetEffTH2Colz = [&](TH2* h,
                                             const string& outPath,
                                             const string& titleText,
                                             const string& zTitle,
                                             bool useLogZ)->void
                {
                  if (!h) return;

                  TCanvas c(TString::Format("c_%s_%s", h->GetName(), rKey.c_str()).Data(), "c_jetEffQA_2D", 950, 800);
                  ApplyCanvasMargins2D(c);
                  c.SetLeftMargin(0.16);
                  c.SetRightMargin(0.16);
                  c.SetBottomMargin(0.14);
                  c.SetTopMargin(0.08);
                  if (useLogZ) c.SetLogz();

                  TH2* hDraw = CloneTH2(h, TString::Format("%s_drawClone", h->GetName()).Data());
                  if (!hDraw) return;
                  hDraw->SetDirectory(nullptr);
                  EnsureSumw2(hDraw);
                  hDraw->SetTitle("");
                  hDraw->GetXaxis()->SetTitle("p_{T}^{#gamma, truth} [GeV]");
                  hDraw->GetYaxis()->SetTitle("x_{J}^{truth}");
                  hDraw->GetZaxis()->SetTitle(zTitle.c_str());
                  hDraw->GetXaxis()->SetTitleOffset(1.10);
                  hDraw->GetYaxis()->SetTitleOffset(1.55);
                  hDraw->GetZaxis()->SetTitleOffset(1.35);

                  if (useLogZ)
                  {
                    double minPos = 1e99;
                    for (int ix = 1; ix <= hDraw->GetNbinsX(); ++ix)
                    {
                      for (int iy = 1; iy <= hDraw->GetNbinsY(); ++iy)
                      {
                        const double z = hDraw->GetBinContent(ix, iy);
                        if (std::isfinite(z) && z > 0.0 && z < minPos) minPos = z;
                      }
                    }
                    if (minPos < 1e98) hDraw->SetMinimum(std::max(0.5 * minPos, 1e-6));
                  }

                  hDraw->Draw("colz");

                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(22);
                  tx.SetTextSize(0.040);
                  tx.DrawLatex(0.50, 0.965, titleText.c_str());

                  SaveCanvas(c, outPath);
                  delete hDraw;
                };

                auto FindExactXBinForAnalysisPt = [&](const TAxis* ax,
                                                      const PtBin& b)->int
                {
                  if (!ax) return -1;

                  for (int ix = 1; ix <= ax->GetNbins(); ++ix)
                  {
                    const double lo = ax->GetBinLowEdge(ix);
                    const double hi = ax->GetBinUpEdge(ix);
                    if (std::fabs(lo - (double)b.lo) < 1e-9 && std::fabs(hi - (double)b.hi) < 1e-9) return ix;
                  }

                  const double cen = 0.5 * (b.lo + b.hi);
                  const int ix = ax->FindBin(cen);
                  if (ix < 1 || ix > ax->GetNbins()) return -1;
                  return ix;
                };

                auto BuildIntegratedVsPtgamma = [&](TH2* h,
                                                    const char* name,
                                                    const char* yTitle)->TH1D*
                {
                  if (!h) return nullptr;

                  vector<double> edges;
                  edges.reserve((std::size_t)nPtAll + 1);
                  for (int i = 0; i < nPtAll; ++i)
                  {
                    if (i == 0) edges.push_back((double)analysisRecoBins[i].lo);
                    edges.push_back((double)analysisRecoBins[i].hi);
                  }

                  TH1D* hOut = new TH1D(name, "", (int)edges.size() - 1, &edges[0]);
                  hOut->SetDirectory(nullptr);
                  hOut->Sumw2();
                  hOut->SetTitle("");
                  hOut->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                  hOut->GetYaxis()->SetTitle(yTitle);

                  for (int i = 0; i < nPtAll; ++i)
                  {
                    const PtBin& b = analysisRecoBins[i];
                    const int ix = FindExactXBinForAnalysisPt(h->GetXaxis(), b);
                    if (ix < 1 || ix > h->GetXaxis()->GetNbins()) continue;

                    double eInt = 0.0;
                    const double val = h->IntegralAndError(
                      ix, ix,
                      1, h->GetYaxis()->GetNbins(),
                      eInt
                    );

                    hOut->SetBinContent(i + 1, val);
                    hOut->SetBinError  (i + 1, eInt);
                  }

                  return hOut;
                };

                auto BuildRatio1D = [&](TH1* hNum,
                                        TH1* hDen,
                                        const char* name,
                                        const char* yTitle)->TH1D*
                {
                  if (!hNum || !hDen) return nullptr;
                  if (hNum->GetNbinsX() != hDen->GetNbinsX()) return nullptr;

                  TH1D* hR = dynamic_cast<TH1D*>(hNum->Clone(name));
                  if (!hR) return nullptr;

                  hR->SetDirectory(nullptr);
                  EnsureSumw2(hR);
                  hR->Reset("ICES");
                  hR->SetTitle("");
                  hR->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                  hR->GetYaxis()->SetTitle(yTitle);

                  for (int ib = 1; ib <= hR->GetNbinsX(); ++ib)
                  {
                    const double num  = hNum->GetBinContent(ib);
                    const double eNum = hNum->GetBinError  (ib);
                    const double den  = hDen->GetBinContent(ib);
                    const double eDen = hDen->GetBinError  (ib);

                    if (!(std::isfinite(num) && std::isfinite(eNum) &&
                          std::isfinite(den) && std::isfinite(eDen) && den > 0.0))
                    {
                      hR->SetBinContent(ib, 0.0);
                      hR->SetBinError  (ib, 0.0);
                      continue;
                    }

                    const double val = num / den;
                    const double var = (eNum * eNum) / (den * den)
                                     + (num * num * eDen * eDen) / (den * den * den * den);
                    const double err = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : 0.0;

                    hR->SetBinContent(ib, val);
                    hR->SetBinError  (ib, err);
                  }

                  return hR;
                };

                auto DrawJetEffIntegrated1D = [&](TH1* h,
                                                  const string& outPath,
                                                  const string& titleText,
                                                  bool drawUnity)->void
                {
                  if (!h) return;

                  TCanvas c(TString::Format("c_%s_%s", h->GetName(), rKey.c_str()).Data(), "c_jetEffQA_1D", 900, 700);
                  ApplyCanvasMargins1D(c);
                  c.SetLeftMargin(0.16);
                  c.SetRightMargin(0.05);
                  c.SetBottomMargin(0.14);
                  c.SetTopMargin(0.08);

                  double maxY = 0.0;
                  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                  {
                    const double v = h->GetBinContent(ib) + h->GetBinError(ib);
                    if (std::isfinite(v) && v > maxY) maxY = v;
                  }

                  TH1F frame("frame_jetEffInt1D", "", 1, 10.0, 35.0);
                  frame.SetMinimum(0.0);
                  frame.SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.05);
                  frame.SetTitle("");
                  frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                  frame.GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
                  frame.GetXaxis()->SetTitleOffset(1.10);
                  frame.GetYaxis()->SetTitleOffset(1.55);
                  frame.Draw("axis");

                  h->SetMarkerStyle(20);
                  h->SetMarkerSize(1.05);
                  h->SetLineWidth(2);
                  h->Draw("E1 same");

                  if (drawUnity)
                  {
                    TLine l1(10.0, 1.0, 35.0, 1.0);
                    l1.SetLineStyle(2);
                    l1.SetLineWidth(2);
                    l1.Draw("same");
                  }

                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(13);
                  tx.SetTextSize(0.034);
                  tx.DrawLatex(0.14, 0.98, titleText.c_str());

                  SaveCanvas(c, outPath);
                };

                auto DrawJetEffIntegratedVsXJ = [&](TH2* hNum,
                                                    TH2* hDen,
                                                    const string& outPath)->void
                {
                  if (!hNum || !hDen) return;

                  const int ny = hDen->GetYaxis()->GetNbins();
                  vector<double> x, ex, y, ey;
                  x.reserve((std::size_t)ny);
                  ex.reserve((std::size_t)ny);
                  y.reserve((std::size_t)ny);
                  ey.reserve((std::size_t)ny);

                  for (int iy = 1; iy <= ny; ++iy)
                  {
                    double eNumInt = 0.0;
                    double eDenInt = 0.0;

                    const double numInt = hNum->IntegralAndError(
                      1, hNum->GetXaxis()->GetNbins(),
                      iy, iy,
                      eNumInt
                    );
                    const double denInt = hDen->IntegralAndError(
                      1, hDen->GetXaxis()->GetNbins(),
                      iy, iy,
                      eDenInt
                    );

                    if (!(denInt > 0.0)) continue;

                    const double xcen = hDen->GetYaxis()->GetBinCenter(iy);
                    const double xerr = 0.5 * hDen->GetYaxis()->GetBinWidth(iy);
                    const double eff  = numInt / denInt;
                    const double var  = (eNumInt * eNumInt) / (denInt * denInt)
                                      + (numInt * numInt * eDenInt * eDenInt) / (denInt * denInt * denInt * denInt);
                    const double err  = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : 0.0;

                    x.push_back(xcen);
                    ex.push_back(xerr);
                    y.push_back(eff);
                    ey.push_back(err);
                  }

                  if (x.empty()) return;

                  TCanvas c(TString::Format("c_jetEffVsXJ_%s", rKey.c_str()).Data(), "c_jetEffVsXJ", 900, 700);
                  ApplyCanvasMargins1D(c);
                  c.SetLeftMargin(0.16);
                  c.SetRightMargin(0.05);
                  c.SetBottomMargin(0.14);
                  c.SetTopMargin(0.08);

                  TH1F frame("frame_jetEffVsXJ", "", 1, 0.0, 2.0);
                  frame.SetMinimum(0.0);
                  frame.SetMaximum(1.05);
                  frame.SetTitle("");
                  frame.GetXaxis()->SetTitle("x_{J}^{truth}");
                  frame.GetYaxis()->SetTitle("Integrated jet efficiency");
                  frame.GetXaxis()->SetTitleOffset(1.10);
                  frame.GetYaxis()->SetTitleOffset(1.55);
                  frame.Draw("axis");

                  TGraphErrors g(
                    (int)x.size(),
                    &x[0], &y[0],
                    &ex[0], &ey[0]
                  );
                  g.SetMarkerStyle(20);
                  g.SetMarkerSize(1.05);
                  g.SetMarkerColor(kBlue + 1);
                  g.SetLineColor(kBlue + 1);
                  g.SetLineWidth(2);
                  g.Draw("P same");

                  TLine l1(0.0, 1.0, 2.0, 1.0);
                  l1.SetLineStyle(2);
                  l1.SetLineWidth(2);
                  l1.Draw("same");

                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(13);
                  tx.SetTextSize(0.034);
                  tx.DrawLatex(
                    0.14, 0.98,
                    TString::Format("Jet Efficiency vs x_{J}^{truth} (R = %.1f), Run24pp, Photon 4 GeV + MBD NS #geq 1", R).Data()
                  );
                  tx.SetTextSize(0.035);
                  tx.DrawLatex(0.14, 0.30, "Integrated over truth p_{T}^{#gamma} unfolding bins");

                  SaveCanvas(c, outPath);
                };

                auto MakeJetEffSliceHist = [&](TH2* hNum,
                                               TH2* hDen,
                                               int ix,
                                               const char* name)->TH1D*
                {
                  if (!hNum || !hDen) return nullptr;
                  if (ix < 1 || ix > hDen->GetXaxis()->GetNbins()) return nullptr;

                  TH1D* hSlice = hDen->ProjectionY(name, ix, ix, "e");
                  if (!hSlice) return nullptr;

                  hSlice->SetDirectory(nullptr);
                  EnsureSumw2(hSlice);
                  hSlice->Reset("ICES");

                  const int ny = hDen->GetYaxis()->GetNbins();
                  for (int iy = 1; iy <= ny; ++iy)
                  {
                    const double num  = hNum->GetBinContent(ix, iy);
                    const double eNum = hNum->GetBinError(ix, iy);
                    const double den  = hDen->GetBinContent(ix, iy);
                    const double eDen = hDen->GetBinError(ix, iy);

                    if (!(std::isfinite(num) && std::isfinite(eNum) &&
                          std::isfinite(den) && std::isfinite(eDen) && den > 0.0))
                    {
                      hSlice->SetBinContent(iy, 0.0);
                      hSlice->SetBinError  (iy, 0.0);
                      continue;
                    }

                    const double val = num / den;
                    const double var = (eNum * eNum) / (den * den)
                                     + (num * num * eDen * eDen) / (den * den * den * den);
                    const double err = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : 0.0;

                    hSlice->SetBinContent(iy, val);
                    hSlice->SetBinError  (iy, err);
                  }

                  hSlice->SetTitle("");
                  hSlice->GetXaxis()->SetTitle("x_{J}^{truth}");
                  hSlice->GetYaxis()->SetTitle("Jet efficiency");
                  hSlice->GetXaxis()->SetTitleOffset(1.10);
                  hSlice->GetYaxis()->SetTitleOffset(1.60);
                  hSlice->SetMarkerStyle(20);
                  hSlice->SetMarkerSize(0.85);
                  hSlice->SetLineWidth(2);
                  return hSlice;
                };

                auto DrawJetEffSliceTable = [&](TH2* hNum,
                                                TH2* hDen,
                                                const string& outPath)->void
                {
                  if (!hNum || !hDen) return;

                  TCanvas c(
                    TString::Format("c_tbl_jetEffSlices_%s", rKey.c_str()).Data(),
                    "c_tbl_jetEffSlices", 1800, 1300
                  );
                  c.Divide(nPtCols, nPtRows, 0.001, 0.001);

                  for (int ipad = 0; ipad < nPtPads; ++ipad)
                  {
                    const int i = ipad;

                    c.cd(ipad + 1);
                    gPad->SetLeftMargin(0.16);
                    gPad->SetRightMargin(0.05);
                    gPad->SetBottomMargin(0.14);
                    gPad->SetTopMargin(0.08);

                    TH1F frame(TString::Format("frame_jetEffSlice_%d", ipad).Data(), "", 1, 0.0, 2.0);
                    frame.SetMinimum(0.0);
                    frame.SetMaximum(1.05);
                    frame.SetTitle("");
                    frame.GetXaxis()->SetTitle("x_{J}^{truth}");
                    frame.GetYaxis()->SetTitle("Jet efficiency");
                    frame.GetXaxis()->SetTitleOffset(1.10);
                    frame.GetYaxis()->SetTitleOffset(1.60);
                    frame.Draw("axis");

                    if (i < 0 || i >= nPtAll)
                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextSize(0.050);
                      tx.DrawLatex(0.16, 0.50, "MISSING");
                      continue;
                    }

                    const PtBin& b = analysisRecoBins[i];
                    const int ix = FindExactXBinForAnalysisPt(hDen->GetXaxis(), b);

                    TH1D* hSlice = MakeJetEffSliceHist(
                      hNum,
                      hDen,
                      ix,
                      TString::Format("h_jetEffSlice_%s_pTbin%d", rKey.c_str(), i + 1).Data()
                    );

                    bool haveSliceContent = false;
                    if (hSlice)
                    {
                      for (int ib = 1; ib <= hSlice->GetNbinsX(); ++ib)
                      {
                        const double y  = hSlice->GetBinContent(ib);
                        const double ey = hSlice->GetBinError(ib);
                        if (std::isfinite(y) && std::isfinite(ey) && (y != 0.0 || ey != 0.0))
                        {
                          haveSliceContent = true;
                          break;
                        }
                      }
                    }

                    if (hSlice && haveSliceContent)
                    {
                      hSlice->GetXaxis()->SetRangeUser(0.0, 2.0);
                      hSlice->SetMinimum(0.0);
                      hSlice->SetMaximum(1.05);
                      hSlice->Draw("E1 same");

                      TLine l1(0.0, 1.0, 2.0, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");
                    }
                    else
                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextSize(0.045);
                      tx.DrawLatex(0.16, 0.50, "NO VALID ENTRIES");
                    }

                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(22);
                      tx.SetTextSize(0.040);
                      tx.DrawLatex(0.52, 0.955,
                                   TString::Format("Jet efficiency vs x_{J}^{truth}, p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                                   b.lo, b.hi, R).Data());
                    }

                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(31);
                      tx.SetTextSize(0.04);
                      const double xR = 0.93;
                      tx.DrawLatex(xR, 0.67, "truth-space efficiency slice");
                      tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                      tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                      tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                    }

                    if (hSlice) delete hSlice;
                  }

                  SaveCanvas(c, outPath);
                };

                auto BuildTH2Ratio = [&](TH2* hNum,
                                         TH2* hDen,
                                         const string& name)->TH2*
                {
                  if (!hNum || !hDen) return nullptr;

                  TH2* hR = CloneTH2(hNum, name);
                  if (!hR) return nullptr;

                  hR->SetDirectory(nullptr);
                  EnsureSumw2(hR);
                  hR->Reset("ICES");

                  for (int ix = 0; ix <= hR->GetNbinsX() + 1; ++ix)
                  {
                    for (int iy = 0; iy <= hR->GetNbinsY() + 1; ++iy)
                    {
                      if (ix == 0 || ix == hR->GetNbinsX() + 1 || iy == 0 || iy == hR->GetNbinsY() + 1)
                      {
                        hR->SetBinContent(ix, iy, 0.0);
                        hR->SetBinError  (ix, iy, 0.0);
                        continue;
                      }

                      const double num  = hNum->GetBinContent(ix, iy);
                      const double eNum = hNum->GetBinError  (ix, iy);
                      const double den  = hDen->GetBinContent(ix, iy);
                      const double eDen = hDen->GetBinError  (ix, iy);

                      if (!(std::isfinite(num) && std::isfinite(eNum) &&
                            std::isfinite(den) && std::isfinite(eDen) && den > 0.0))
                      {
                        hR->SetBinContent(ix, iy, 0.0);
                        hR->SetBinError  (ix, iy, 0.0);
                        continue;
                      }

                      const double val = num / den;
                      const double var = (eNum * eNum) / (den * den)
                                       + (num * num * eDen * eDen) / (den * den * den * den);
                      const double err = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : 0.0;

                      hR->SetBinContent(ix, iy, val);
                      hR->SetBinError  (ix, iy, err);
                    }
                  }

                  return hR;
                };

                h2JetEff = CloneTH2(
                  h2JetEffNum_in,
                  TString::Format("h2JetEff_%s", rKey.c_str()).Data()
                );

                if (h2JetEff)
                {
                  h2JetEff->SetDirectory(nullptr);
                  EnsureSumw2(h2JetEff);
                  h2JetEff->Divide(h2JetEffDen_in);

                  vector<double> xPt, exPt, yEff, eyEff;
                  xPt.reserve((std::size_t)nPtAll);
                  exPt.reserve((std::size_t)nPtAll);
                  yEff.reserve((std::size_t)nPtAll);
                  eyEff.reserve((std::size_t)nPtAll);

                  for (int i = 0; i < nPtAll; ++i)
                  {
                    const PtBin& b = analysisRecoBins[i];
                    const double ex  = 0.5 * (b.hi - b.lo);

                    const int ixEff = FindExactXBinForAnalysisPt(h2JetEffDen_in->GetXaxis(), b);
                    if (ixEff < 1 || ixEff > h2JetEffDen_in->GetXaxis()->GetNbins()) continue;

                    double eNumInt = 0.0;
                    double eDenInt = 0.0;

                    const double numInt = h2JetEffNum_in->IntegralAndError(
                      ixEff, ixEff,
                      1, h2JetEffNum_in->GetYaxis()->GetNbins(),
                      eNumInt
                    );
                    const double denInt = h2JetEffDen_in->IntegralAndError(
                      ixEff, ixEff,
                      1, h2JetEffDen_in->GetYaxis()->GetNbins(),
                      eDenInt
                    );

                    if (!(denInt > 0.0)) continue;

                    const double eff = numInt / denInt;
                    const double var = (eNumInt * eNumInt) / (denInt * denInt)
                                     + (numInt * numInt * eDenInt * eDenInt) / (denInt * denInt * denInt * denInt);
                    const double err = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : 0.0;

                    xPt.push_back(0.5 * (b.lo + b.hi));
                    exPt.push_back(ex);
                    yEff.push_back(eff);
                    eyEff.push_back(err);
                  }

                  if (!xPt.empty())
                  {
                    TCanvas cJetEffPt(
                      TString::Format("c_jetEffPt_%s", rKey.c_str()).Data(),
                      "c_jetEffPt", 900, 700
                    );
                    ApplyCanvasMargins1D(cJetEffPt);
                    cJetEffPt.SetLeftMargin(0.16);
                    cJetEffPt.SetRightMargin(0.05);
                    cJetEffPt.SetBottomMargin(0.14);
                    cJetEffPt.SetTopMargin(0.08);

                    TH1F frame("frame_jetEffPt","", 1, 10.0, 35.0);
                    frame.SetMinimum(0.0);
                    frame.SetMaximum(1.05);
                    frame.SetTitle("");
                    frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    frame.GetYaxis()->SetTitle("Integrated jet efficiency");
                    frame.GetXaxis()->SetTitleOffset(1.10);
                    frame.GetYaxis()->SetTitleOffset(1.55);
                    frame.Draw("axis");

                    TGraphErrors gJetEff(
                      (int)xPt.size(),
                      &xPt[0], &yEff[0],
                      &exPt[0], &eyEff[0]
                    );
                    gJetEff.SetMarkerStyle(20);
                    gJetEff.SetMarkerSize(1.05);
                    gJetEff.SetMarkerColor(kBlue + 1);
                    gJetEff.SetLineColor(kBlue + 1);
                    gJetEff.SetLineWidth(2);
                    gJetEff.Draw("P same");

                    TLine l1(10.0, 1.0, 35.0, 1.0);
                    l1.SetLineStyle(2);
                    l1.SetLineWidth(2);
                    l1.Draw("same");

                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(13);
                      tx.SetTextSize(0.034);
                      tx.DrawLatex(
                        0.14, 0.98,
                        TString::Format("Jet Efficiency vs p_{T}^{#gamma} (R = %.1f), Run24pp, Photon 4 GeV + MBD NS #geq 1", R).Data()
                      );
                    }

                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(13);
                      tx.SetTextSize(0.035);
                      tx.DrawLatex(0.14, 0.30, "Integrated over truth x_{J} unfolding bins");
                    }

                    SaveCanvas(cJetEffPt, JoinPath(withAndWithoutJetEffOut, "jetEfficiency_integrated_vs_pTgamma.png"));
                  }

                  DrawJetEffTH2Colz(
                    h2JetEffDen_in,
                    JoinPath(jetEffQAOut, "jetEff_truthPhaseSpace_den_pTgamma_xJ_colz.png"),
                    TString::Format("Jet-efficiency denominator, R = %.1f", R).Data(),
                    "Truth recoil-jet counts",
                    true
                  );

                  DrawJetEffTH2Colz(
                    h2JetEffNum_in,
                    JoinPath(jetEffQAOut, "jetEff_truthPhaseSpace_num_pTgamma_xJ_colz.png"),
                    TString::Format("Jet-efficiency numerator, R = %.1f", R).Data(),
                    "Matched truth recoil-jet counts",
                    true
                  );

                  DrawJetEffTH2Colz(
                    h2JetEff,
                    JoinPath(jetEffQAOut, "jetEff_truthPhaseSpace_eff_pTgamma_xJ_colz.png"),
                    TString::Format("Jet efficiency map, R = %.1f", R).Data(),
                    "Jet efficiency",
                    false
                  );

                  DrawJetEffIntegratedVsXJ(
                    h2JetEffNum_in,
                    h2JetEffDen_in,
                    JoinPath(jetEffQAOut, "jetEfficiency_integrated_vs_xJ.png")
                  );

                  DrawJetEffSliceTable(
                    h2JetEffNum_in,
                    h2JetEffDen_in,
                    JoinPath(jetEffQAOut, "table3x3_jetEfficiency_vs_xJ_truthSlices.png")
                  );

                  TH1D* hTruthJetIntegrated = BuildIntegratedVsPtgamma(
                    h2JetEffDen_in,
                    TString::Format("hTruthJetIntegrated_%s", rKey.c_str()).Data(),
                    "Integrated truth recoil-jet counts"
                  );
                  TH1D* hRecoJetMCIntegrated = BuildIntegratedVsPtgamma(
                    h2RecoSim,
                    TString::Format("hRecoJetMCIntegrated_%s", rKey.c_str()).Data(),
                    "Integrated reco-jet counts (MC)"
                  );
                  TH1D* hRecoJetDataIntegrated = BuildIntegratedVsPtgamma(
                    h2RecoData,
                    TString::Format("hRecoJetDataIntegrated_%s", rKey.c_str()).Data(),
                    "Integrated reco-jet counts (DATA)"
                  );

                  TH1D* hTruthOverRecoMC = BuildRatio1D(
                    hTruthJetIntegrated,
                    hRecoJetMCIntegrated,
                    TString::Format("hTruthOverRecoMC_%s", rKey.c_str()).Data(),
                    "Truth jet / reco jet (MC)"
                  );
                  TH1D* hTruthOverData = BuildRatio1D(
                    hTruthJetIntegrated,
                    hRecoJetDataIntegrated,
                    TString::Format("hTruthOverData_%s", rKey.c_str()).Data(),
                    "Truth jet / reco jet (DATA)"
                  );
                  TH1D* hRecoMCOverData = BuildRatio1D(
                    hRecoJetMCIntegrated,
                    hRecoJetDataIntegrated,
                    TString::Format("hRecoMCOverData_%s", rKey.c_str()).Data(),
                    "Reco jet (MC) / reco jet (DATA)"
                  );

                  if (hTruthOverRecoMC)
                  {
                    hTruthOverRecoMC->SetMarkerColor(kBlack);
                    hTruthOverRecoMC->SetLineColor(kBlack);
                    DrawJetEffIntegrated1D(
                      hTruthOverRecoMC,
                      JoinPath(jetEffQAOut, "jetTruthOverReco_MC_integrated_vs_pTgamma.png"),
                      TString::Format("Truth / reco jet ratio (MC), R = %.1f", R).Data(),
                      true
                    );
                  }

                  if (hTruthOverData)
                  {
                    hTruthOverData->SetMarkerColor(kRed + 1);
                    hTruthOverData->SetLineColor(kRed + 1);
                    DrawJetEffIntegrated1D(
                      hTruthOverData,
                      JoinPath(jetEffQAOut, "jetTruthOverReco_DATA_integrated_vs_pTgamma.png"),
                      TString::Format("Truth / reco jet ratio (truth vs DATA reco), R = %.1f", R).Data(),
                      true
                    );
                  }

                  if (hRecoMCOverData)
                  {
                    hRecoMCOverData->SetMarkerColor(kBlue + 1);
                    hRecoMCOverData->SetLineColor(kBlue + 1);
                    DrawJetEffIntegrated1D(
                      hRecoMCOverData,
                      JoinPath(jetEffQAOut, "jetRecoMCOverRecoDATA_integrated_vs_pTgamma.png"),
                      TString::Format("Reco jet ratio (MC / DATA), R = %.1f", R).Data(),
                      true
                    );
                  }

                  if (h2UnfoldTruth)
                  {
                    bool sameJetEffBinning = true;

                    if (h2JetEff->GetNbinsX() != h2UnfoldTruth->GetNbinsX() ||
                        h2JetEff->GetNbinsY() != h2UnfoldTruth->GetNbinsY())
                    {
                      sameJetEffBinning = false;
                    }

                    if (sameJetEffBinning)
                    {
                      for (int ix = 1; ix <= h2UnfoldTruth->GetNbinsX() + 1; ++ix)
                      {
                        const double e1 = h2JetEff->GetXaxis()->GetBinUpEdge(ix);
                        const double e2 = h2UnfoldTruth->GetXaxis()->GetBinUpEdge(ix);
                        if (std::fabs(e1 - e2) > 1e-9)
                        {
                          sameJetEffBinning = false;
                          break;
                        }
                      }
                    }

                    if (sameJetEffBinning)
                    {
                      for (int iy = 1; iy <= h2UnfoldTruth->GetNbinsY() + 1; ++iy)
                      {
                        const double e1 = h2JetEff->GetYaxis()->GetBinUpEdge(iy);
                        const double e2 = h2UnfoldTruth->GetYaxis()->GetBinUpEdge(iy);
                        if (std::fabs(e1 - e2) > 1e-9)
                        {
                          sameJetEffBinning = false;
                          break;
                        }
                      }
                    }

                    if (!sameJetEffBinning)
                    {
                      cout << ANSI_BOLD_YEL
                           << "[WARN] Jet-efficiency binning mismatch for " << rKey
                           << ". Skipping 2D jet-efficiency correction for this radius."
                           << ANSI_RESET << "\n";
                    }
                    else
                    {
                      h2UnfoldTruth_jetEffCorr = CloneTH2(
                        h2UnfoldTruth,
                        TString::Format("h2_unfoldedTruth_pTgamma_xJ_incl_jetEffCorr_%s", rKey.c_str()).Data()
                      );

                      if (h2UnfoldTruth_jetEffCorr)
                      {
                        h2UnfoldTruth_jetEffCorr->SetDirectory(nullptr);
                        EnsureSumw2(h2UnfoldTruth_jetEffCorr);
                        h2UnfoldTruth_jetEffCorr->Reset("ICES");

                        const int nxJC = h2UnfoldTruth->GetNbinsX();
                        const int nyJC = h2UnfoldTruth->GetNbinsY();

                        for (int ix = 0; ix <= nxJC + 1; ++ix)
                        {
                          for (int iy = 0; iy <= nyJC + 1; ++iy)
                          {
                            if (ix == 0 || ix == nxJC + 1 || iy == 0 || iy == nyJC + 1)
                            {
                              h2UnfoldTruth_jetEffCorr->SetBinContent(ix, iy, 0.0);
                              h2UnfoldTruth_jetEffCorr->SetBinError  (ix, iy, 0.0);
                              continue;
                            }

                            const double num  = h2UnfoldTruth->GetBinContent(ix, iy);
                            const double eNum = h2UnfoldTruth->GetBinError  (ix, iy);
                            const double eff  = h2JetEff->GetBinContent(ix, iy);
                            const double eEff = h2JetEff->GetBinError  (ix, iy);

                            if (!(std::isfinite(num) && std::isfinite(eNum) &&
                                  std::isfinite(eff) && std::isfinite(eEff) && eff > 0.0))
                            {
                              h2UnfoldTruth_jetEffCorr->SetBinContent(ix, iy, 0.0);
                              h2UnfoldTruth_jetEffCorr->SetBinError  (ix, iy, 0.0);
                              continue;
                            }

                            const double corr = num / eff;
                            const double var  = (eNum * eNum) / (eff * eff)
                                              + (num * num * eEff * eEff) / (eff * eff * eff * eff);
                            const double err  = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : 0.0;

                            h2UnfoldTruth_jetEffCorr->SetBinContent(ix, iy, corr);
                            h2UnfoldTruth_jetEffCorr->SetBinError  (ix, iy, err);
                          }
                        }

                        DrawJetEffTH2Colz(
                          h2UnfoldTruth,
                          JoinPath(jetEffQAOut, "unfoldedTruth_withoutJetEffCorr_pTgamma_xJ_colz.png"),
                          TString::Format("Unfolded truth without jet-eff. correction, R = %.1f", R).Data(),
                          "Unfolded yield",
                          true
                        );

                        DrawJetEffTH2Colz(
                          h2UnfoldTruth_jetEffCorr,
                          JoinPath(jetEffQAOut, "unfoldedTruth_withJetEffCorr_pTgamma_xJ_colz.png"),
                          TString::Format("Unfolded truth with jet-eff. correction, R = %.1f", R).Data(),
                          "Corrected unfolded yield",
                          true
                        );

                        TH2* h2CorrOverUncorr = BuildTH2Ratio(
                          h2UnfoldTruth_jetEffCorr,
                          h2UnfoldTruth,
                          TString::Format("h2_corrOverUncorr_%s", rKey.c_str()).Data()
                        );

                        if (h2CorrOverUncorr)
                        {
                          DrawJetEffTH2Colz(
                            h2CorrOverUncorr,
                            JoinPath(jetEffQAOut, "unfoldedTruth_ratio_withOverWithoutJetEffCorr_pTgamma_xJ_colz.png"),
                            TString::Format("Effect of jet-eff. correction, R = %.1f", R).Data(),
                            "With jet-eff. corr. / without",
                            false
                          );
                          delete h2CorrOverUncorr;
                        }
                      }
                    }
                  }

                  if (hTruthOverRecoMC) delete hTruthOverRecoMC;
                  if (hTruthOverData) delete hTruthOverData;
                  if (hRecoMCOverData) delete hRecoMCOverData;
                  if (hTruthJetIntegrated) delete hTruthJetIntegrated;
                  if (hRecoJetMCIntegrated) delete hRecoJetMCIntegrated;
                  if (hRecoJetDataIntegrated) delete hRecoJetDataIntegrated;
                }
             }
          }

          auto BuildRatioHist = [&](TH1* hNum, TH1* hDen, const char* newName)->TH1*
        {
            if (!hNum || !hDen) return nullptr;

            TH1* hR = CloneTH1(hNum, newName);
            if (!hR) return nullptr;

            hR->SetDirectory(nullptr);
            EnsureSumw2(hR);
            hR->Reset("ICES");

            const int nb = hR->GetNbinsX();
            for (int ib = 0; ib <= nb + 1; ++ib)
            {
              if (ib == 0 || ib == nb + 1)
              {
                hR->SetBinContent(ib, 0.0);
                hR->SetBinError  (ib, 0.0);
                continue;
              }

              const double num  = hNum->GetBinContent(ib);
              const double eNum = hNum->GetBinError  (ib);
              const double den  = hDen->GetBinContent(ib);
              const double eDen = hDen->GetBinError  (ib);

              if (!(std::isfinite(num) && std::isfinite(eNum) && std::isfinite(den) && std::isfinite(eDen)) || den <= 0.0)
              {
                hR->SetBinContent(ib, 0.0);
                hR->SetBinError  (ib, 0.0);
                continue;
              }

              const double val = num / den;
              const double var = (eNum * eNum) / (den * den)
                               + (num * num * eDen * eDen) / (den * den * den * den);
              const double err = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : 0.0;

              hR->SetBinContent(ib, val);
              hR->SetBinError  (ib, err);
            }

            hR->SetTitle("");
            hR->GetXaxis()->SetTitle("x_{J}");
            return hR;
        };
        // ----------------------------------------------------------------------
        // Closure test (SIM): truth -> smear(reco) -> unfold -> compare back to truth
        //   Summary vs pT: (Integral of unfolded xJ)/(Integral of truth xJ)  ~ 1
        //   Output: <rOut>/closure_unfoldedOverTruth_integral_vs_pTgamma.png
        // ----------------------------------------------------------------------
          {
              TH1D* hMeasSimGlob_closure = hRsp_measXtruth->ProjectionX(
                TString::Format("hMeasSimGlob_forClosure_%s", rKey.c_str()).Data(),
                0, hRsp_measXtruth->GetNbinsY() + 1, "e"
              );
              TH1D* hTruthSimGlob_closure = hRsp_measXtruth->ProjectionY(
                TString::Format("hTruthSimGlob_forClosure_%s", rKey.c_str()).Data(),
                0, hRsp_measXtruth->GetNbinsX() + 1, "e"
              );
              if (hMeasSimGlob_closure)  { hMeasSimGlob_closure->SetDirectory(nullptr);  EnsureSumw2(hMeasSimGlob_closure); }
              if (hTruthSimGlob_closure) { hTruthSimGlob_closure->SetDirectory(nullptr); EnsureSumw2(hTruthSimGlob_closure); }

              RooUnfoldResponse respXJ_closure(
                hMeasSimGlob_closure, hTruthSimGlob_closure, hRsp_measXtruth,
                TString::Format("respXJ_closure_%s", rKey.c_str()).Data(),
                TString::Format("respXJ_closure_%s", rKey.c_str()).Data()
              );

              RooUnfoldBayes unfoldXJ_closure(
                &respXJ_closure,
                (hMeasSimGlob_closure ? hMeasSimGlob_closure : hMeasSimGlob),
                kBayesIterXJ
              );
              unfoldXJ_closure.SetVerbose(0);
              unfoldXJ_closure.SetNToys(kNToysXJScan);

              TH1* hUnfoldTruthGlob_closure = nullptr;
              if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
              hUnfoldTruthGlob_closure = unfoldXJ_closure.Hreco(RooUnfold::kCovToy);
              if (gSystem) gSystem->RedirectOutput(0);
              if (hUnfoldTruthGlob_closure) hUnfoldTruthGlob_closure->SetDirectory(nullptr);

              TH2* h2UnfoldTruth_closure = nullptr;
              if (hUnfoldTruthGlob_closure)
              {
                h2UnfoldTruth_closure = UnflattenGlobalToTH2(
                  hUnfoldTruthGlob_closure,
                  h2TruthSim,
                  TString::Format("h2_unfoldedTruth_closure_pTgamma_xJ_%s", rKey.c_str()).Data()
                );
              }

              TH2* h2Truth_closure = nullptr;
              if (hTruthSimGlob_closure)
              {
                h2Truth_closure = UnflattenGlobalToTH2(
                  hTruthSimGlob_closure,
                  h2TruthSim,
                  TString::Format("h2_truthClosure_pTgamma_xJ_%s", rKey.c_str()).Data()
                );
              }

              vector<double> xPt, exPt, yRat, eyRat;

              const auto& analysisRecoBins = UnfoldAnalysisRecoPtBins();
              const int nPtClosure = (int)analysisRecoBins.size();
              for (int i = 0; i < nPtClosure; ++i)
              {
                const PtBin& b = analysisRecoBins[i];
                const double cen = 0.5 * (b.lo + b.hi);
                const double ex  = 0.5 * (b.hi - b.lo);

                if (!h2UnfoldTruth_closure || !h2Truth_closure) continue;

                const int ixTruth = h2Truth_closure->GetXaxis()->FindBin(cen);
                if (ixTruth < 1 || ixTruth > h2Truth_closure->GetXaxis()->GetNbins()) continue;

                TH1D* hXJ_truth = h2Truth_closure->ProjectionY(
                  TString::Format("h_xJ_truthClosure_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
                  ixTruth, ixTruth, "e"
                );
                TH1D* hXJ_unf = h2UnfoldTruth_closure->ProjectionY(
                  TString::Format("h_xJ_unfClosure_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
                  ixTruth, ixTruth, "e"
                );

                if (!hXJ_truth || !hXJ_unf)
                {
                  if (hXJ_truth) delete hXJ_truth;
                  if (hXJ_unf)   delete hXJ_unf;
                  continue;
                }

                hXJ_truth->SetDirectory(nullptr);
                hXJ_unf->SetDirectory(nullptr);
                EnsureSumw2(hXJ_truth);
                EnsureSumw2(hXJ_unf);

                double eTruth = 0.0;
                double eUnf   = 0.0;

                const double ITruth = hXJ_truth->IntegralAndError(1, hXJ_truth->GetNbinsX(), eTruth, "width");
                const double IUnf   = hXJ_unf  ->IntegralAndError(1, hXJ_unf  ->GetNbinsX(), eUnf,   "width");

                if (ITruth > 0.0)
                {
                  const double r  = IUnf / ITruth;

                  double relUnf = 0.0;
                  if (IUnf > 0.0) relUnf = eUnf / IUnf;

                  const double relTruth = eTruth / ITruth;
                  const double er = std::fabs(r) * std::sqrt(relUnf*relUnf + relTruth*relTruth);

                  xPt.push_back(cen);
                  exPt.push_back(ex);
                  yRat.push_back(r);
                  eyRat.push_back(er);
                }

                delete hXJ_truth;
                delete hXJ_unf;
              }

              if (!xPt.empty())
              {
                  const double ymin = 0.95;
                  const double ymax = 1.05;

                  // Save the r04 closure y-range to be reused by photon closure plots
                  if (rKey == "r04")
                  {
                      gYmin_r04_closure = ymin;
                      gYmax_r04_closure = ymax;
                  }

                  TCanvas c(TString::Format("c_closure_vs_pt_%s", rKey.c_str()).Data(), "c_closure_vs_pt", 900, 700);
                  ApplyCanvasMargins1D(c);

                  TGraphErrors g((int)xPt.size(), &xPt[0], &yRat[0], &exPt[0], &eyRat[0]);
                  g.SetLineWidth(2);
                  g.SetMarkerStyle(20);
                  g.SetTitle("");
                  g.Draw("AP");

                  g.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                  g.GetXaxis()->SetTitleOffset(1.15);
                  g.GetYaxis()->SetTitle("Closure: Unfolded MC / Truth MC");
                  g.GetYaxis()->SetRangeUser(ymin, ymax);

                  TLine l1(g.GetXaxis()->GetXmin(), 1.0, g.GetXaxis()->GetXmax(), 1.0);
                  l1.SetLineStyle(2);
                  l1.SetLineWidth(2);
                  l1.Draw("same");

                  // Centered title (replaces the old TLatex header lines)
                  {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(22);
                    tx.SetTextSize(0.040);
                    tx.DrawLatex(0.50, 0.965,
                                 TString::Format("Closure test, unfold(reco) #rightarrow truth, R = %.1f", R).Data());
                  }

                  // Top-left Bayes annotation (moved up/left)
                  {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(13);
                    tx.SetTextSize(0.038);
                    tx.DrawLatex(0.20, 0.90, TString::Format("Bayes it=%d (xJ)", kBayesIterXJ).Data());
                    tx.DrawLatex(0.20, 0.855, "2D (p_{T}^{#gamma}, x_{J}) unfolding");
                  }


                  SaveCanvas(c, JoinPath(closurePlotsOut, "closure_unfoldedOverTruth_integral_vs_pTgamma.png"));
              }

              if (h2Truth_closure) delete h2Truth_closure;
              if (h2UnfoldTruth_closure) delete h2UnfoldTruth_closure;
              if (hUnfoldTruthGlob_closure) delete hUnfoldTruthGlob_closure;
              if (hTruthSimGlob_closure) delete hTruthSimGlob_closure;
              if (hMeasSimGlob_closure) delete hMeasSimGlob_closure;
          }

        // ----------------------------------------------------------------------
        // Half-closure test (SIM, single-file infrastructure):
        //   Split the SIM response matrix into two disjoint halves (A/B) at the
        //   response-matrix bin level (Binomial partition, 50/50).
        //
        //   Train response on half A, unfold half B reco, compare to half B truth.
        //
        //   Summary vs pT: (Integral of unfolded_B xJ)/(Integral of truth_B xJ)  ~ 1
        //
        //   Output: <rOut>/halfClosure_unfoldedOverTruth_integral_vs_pTgamma.png
        // ----------------------------------------------------------------------
        {
            if (hRsp_measXtruth && hMeasSimGlob && hTruthSimGlob && h2TruthSim)
            {
              unsigned int seed = 1337u;
              for (size_t ic = 0; ic < rKey.size(); ++ic)
              {
                seed = seed * 131u + (unsigned int)(unsigned char)rKey[ic];
              }
              std::mt19937 rng(seed);

              TH2D* hRsp_measXtruth_A = dynamic_cast<TH2D*>(CloneTH2(
                hRsp_measXtruth,
                TString::Format("h2_unfoldResponse_global_recoVsTruth_halfA_%s", rKey.c_str()).Data()
              ));
              TH2D* hRsp_measXtruth_B = dynamic_cast<TH2D*>(CloneTH2(
                hRsp_measXtruth,
                TString::Format("h2_unfoldResponse_global_recoVsTruth_halfB_%s", rKey.c_str()).Data()
              ));

              TH1* hMeasSimGlob_A  = CloneTH1(hMeasSimGlob,
                TString::Format("hMeasSimGlob_halfA_%s", rKey.c_str()).Data());
              TH1* hTruthSimGlob_A = CloneTH1(hTruthSimGlob,
                TString::Format("hTruthSimGlob_halfA_%s", rKey.c_str()).Data());
              TH1* hMeasSimGlob_B  = CloneTH1(hMeasSimGlob,
                TString::Format("hMeasSimGlob_halfB_%s", rKey.c_str()).Data());
              TH1* hTruthSimGlob_B = CloneTH1(hTruthSimGlob,
                TString::Format("hTruthSimGlob_halfB_%s", rKey.c_str()).Data());

              if (hRsp_measXtruth_A) { hRsp_measXtruth_A->Reset("ICES"); hRsp_measXtruth_A->SetDirectory(nullptr); hRsp_measXtruth_A->Sumw2(); }
              if (hRsp_measXtruth_B) { hRsp_measXtruth_B->Reset("ICES"); hRsp_measXtruth_B->SetDirectory(nullptr); hRsp_measXtruth_B->Sumw2(); }

              if (hMeasSimGlob_A)  { hMeasSimGlob_A ->Reset("ICES"); EnsureSumw2(hMeasSimGlob_A); }
              if (hTruthSimGlob_A) { hTruthSimGlob_A->Reset("ICES"); EnsureSumw2(hTruthSimGlob_A); }
              if (hMeasSimGlob_B)  { hMeasSimGlob_B ->Reset("ICES"); EnsureSumw2(hMeasSimGlob_B); }
              if (hTruthSimGlob_B) { hTruthSimGlob_B->Reset("ICES"); EnsureSumw2(hTruthSimGlob_B); }

              if (!hRsp_measXtruth_A || !hRsp_measXtruth_B || !hMeasSimGlob_A || !hTruthSimGlob_A || !hMeasSimGlob_B || !hTruthSimGlob_B)
              {
                cout << ANSI_BOLD_YEL
                     << "[WARN] Half-closure: failed to allocate split response/marginals for " << rKey
                     << ". Skipping half-closure."
                     << ANSI_RESET << "\n";
              }
              else
              {
                vector<double> measA((std::size_t)nGlobReco_rsp  + 2, 0.0);
                vector<double> measB((std::size_t)nGlobReco_rsp  + 2, 0.0);
                vector<double> truthA((std::size_t)nGlobTruth_rsp + 2, 0.0);
                vector<double> truthB((std::size_t)nGlobTruth_rsp + 2, 0.0);

                for (int ixTruth = 0; ixTruth <= nGlobTruth_rsp + 1; ++ixTruth)
                {
                  for (int iyReco = 0; iyReco <= nGlobReco_rsp + 1; ++iyReco)
                  {
                    const double nRaw = hRsp_measXtruth->GetBinContent(iyReco, ixTruth);
                    if (!(nRaw > 0.0))
                    {
                      hRsp_measXtruth_A->SetBinContent(iyReco, ixTruth, 0.0);
                      hRsp_measXtruth_A->SetBinError  (iyReco, ixTruth, 0.0);
                      hRsp_measXtruth_B->SetBinContent(iyReco, ixTruth, 0.0);
                      hRsp_measXtruth_B->SetBinError  (iyReco, ixTruth, 0.0);
                      continue;
                    }

                    const long long N = (long long)std::llround(nRaw);
                    if (N <= 0LL)
                    {
                      hRsp_measXtruth_A->SetBinContent(iyReco, ixTruth, 0.0);
                      hRsp_measXtruth_A->SetBinError  (iyReco, ixTruth, 0.0);
                      hRsp_measXtruth_B->SetBinContent(iyReco, ixTruth, 0.0);
                      hRsp_measXtruth_B->SetBinError  (iyReco, ixTruth, 0.0);
                      continue;
                    }

                    std::binomial_distribution<long long> d(N, 0.5);
                    const long long NA = d(rng);
                    const long long NB = N - NA;

                    hRsp_measXtruth_A->SetBinContent(iyReco, ixTruth, (double)NA);
                    hRsp_measXtruth_A->SetBinError  (iyReco, ixTruth, std::sqrt((double)NA));

                    hRsp_measXtruth_B->SetBinContent(iyReco, ixTruth, (double)NB);
                    hRsp_measXtruth_B->SetBinError  (iyReco, ixTruth, std::sqrt((double)NB));

                    if (iyReco >= 0 && iyReco <= nGlobReco_rsp + 1)
                    {
                      measA[(std::size_t)iyReco] += (double)NA;
                      measB[(std::size_t)iyReco] += (double)NB;
                    }
                    if (ixTruth >= 0 && ixTruth <= nGlobTruth_rsp + 1)
                    {
                      truthA[(std::size_t)ixTruth] += (double)NA;
                      truthB[(std::size_t)ixTruth] += (double)NB;
                    }
                  }
                }

                for (int ib = 0; ib <= nGlobReco_rsp + 1; ++ib)
                {
                  hMeasSimGlob_A->SetBinContent(ib, measA[(std::size_t)ib]);
                  hMeasSimGlob_A->SetBinError  (ib, std::sqrt(measA[(std::size_t)ib]));
                  hMeasSimGlob_B->SetBinContent(ib, measB[(std::size_t)ib]);
                  hMeasSimGlob_B->SetBinError  (ib, std::sqrt(measB[(std::size_t)ib]));
                }

                for (int ib = 0; ib <= nGlobTruth_rsp + 1; ++ib)
                {
                  hTruthSimGlob_A->SetBinContent(ib, truthA[(std::size_t)ib]);
                  hTruthSimGlob_A->SetBinError  (ib, std::sqrt(truthA[(std::size_t)ib]));
                  hTruthSimGlob_B->SetBinContent(ib, truthB[(std::size_t)ib]);
                  hTruthSimGlob_B->SetBinError  (ib, std::sqrt(truthB[(std::size_t)ib]));
                }

                RooUnfoldResponse respXJ_half(hMeasSimGlob_A, hTruthSimGlob_A, hRsp_measXtruth_A,
                                             TString::Format("respXJ_halfA_%s", rKey.c_str()).Data(),
                                             TString::Format("respXJ_halfA_%s", rKey.c_str()).Data());

                TH1* hMeasSimGlob_halfB_meas = CloneTH1(hMeasSimGlob_B,
                  TString::Format("hMeasSimGlob_forHalfClosure_%s", rKey.c_str()).Data());
                if (hMeasSimGlob_halfB_meas) hMeasSimGlob_halfB_meas->SetDirectory(nullptr);

                RooUnfoldBayes unfoldXJ_half(&respXJ_half, (hMeasSimGlob_halfB_meas ? hMeasSimGlob_halfB_meas : hMeasSimGlob_B), kBayesIterXJ);
                  unfoldXJ_half.SetVerbose(0);
                  unfoldXJ_half.SetNToys(kNToysXJScan);

                TH1* hUnfoldTruthGlob_half = nullptr;
                if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
                hUnfoldTruthGlob_half = unfoldXJ_half.Hreco(RooUnfold::kCovToy);
                if (gSystem) gSystem->RedirectOutput(0);
                if (hUnfoldTruthGlob_half) hUnfoldTruthGlob_half->SetDirectory(nullptr);

                TH2* h2UnfoldTruth_half = nullptr;
                if (hUnfoldTruthGlob_half)
                {
                  h2UnfoldTruth_half = UnflattenGlobalToTH2(
                    hUnfoldTruthGlob_half,
                    h2TruthSim,
                    TString::Format("h2_unfoldedTruth_halfClosure_pTgamma_xJ_%s", rKey.c_str()).Data()
                  );
                }

                TH2* h2Truth_halfB = nullptr;
                if (hTruthSimGlob_B)
                {
                  h2Truth_halfB = UnflattenGlobalToTH2(
                    hTruthSimGlob_B,
                    h2TruthSim,
                    TString::Format("h2_truth_halfB_pTgamma_xJ_%s", rKey.c_str()).Data()
                  );
                }

                vector<double> xPt, exPt, yRat, eyRat;

                const auto& analysisRecoBins = UnfoldAnalysisRecoPtBins();
                const int nPtHalf = (int)analysisRecoBins.size();
                for (int i = 0; i < nPtHalf; ++i)
                {
                  const PtBin& b = analysisRecoBins[i];
                  const double cen = 0.5 * (b.lo + b.hi);
                  const double ex  = 0.5 * (b.hi - b.lo);

                  if (!h2UnfoldTruth_half || !h2Truth_halfB) continue;

                  const int ixTruth = h2Truth_halfB->GetXaxis()->FindBin(cen);
                  if (ixTruth < 1 || ixTruth > h2Truth_halfB->GetXaxis()->GetNbins()) continue;

                  TH1D* hXJ_truth = h2Truth_halfB->ProjectionY(
                    TString::Format("h_xJ_truthHalfClosure_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
                    ixTruth, ixTruth, "e"
                  );
                  TH1D* hXJ_unf = h2UnfoldTruth_half->ProjectionY(
                    TString::Format("h_xJ_unfHalfClosure_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
                    ixTruth, ixTruth, "e"
                  );

                  if (!hXJ_truth || !hXJ_unf)
                  {
                    if (hXJ_truth) delete hXJ_truth;
                    if (hXJ_unf)   delete hXJ_unf;
                    continue;
                  }

                  hXJ_truth->SetDirectory(nullptr);
                  hXJ_unf->SetDirectory(nullptr);
                  EnsureSumw2(hXJ_truth);
                  EnsureSumw2(hXJ_unf);

                  double eTruth = 0.0;
                  double eUnf   = 0.0;

                  const double ITruth = hXJ_truth->IntegralAndError(1, hXJ_truth->GetNbinsX(), eTruth, "width");
                  const double IUnf   = hXJ_unf  ->IntegralAndError(1, hXJ_unf  ->GetNbinsX(), eUnf,   "width");

                  if (ITruth > 0.0)
                  {
                    const double r  = IUnf / ITruth;

                    double relUnf = 0.0;
                    if (IUnf > 0.0) relUnf = eUnf / IUnf;

                    const double relTruth = eTruth / ITruth;
                    const double er = std::fabs(r) * std::sqrt(relUnf*relUnf + relTruth*relTruth);

                    xPt.push_back(cen);
                    exPt.push_back(ex);
                    yRat.push_back(r);
                    eyRat.push_back(er);
                  }

                  delete hXJ_truth;
                  delete hXJ_unf;
                }

                if (!xPt.empty())
                {
                  double ymin =  1e99;
                  double ymax = -1e99;
                  for (size_t i = 0; i < yRat.size(); ++i)
                  {
                    ymin = std::min(ymin, yRat[i] - eyRat[i]);
                    ymax = std::max(ymax, yRat[i] + eyRat[i]);
                  }
                    if (!(ymin < 1e98) || !(ymax > -1e98) || ymin >= ymax)
                    {
                      ymin = 0.5;
                      ymax = 1.5;
                    }
                    else
                    {
                      const double pad = 0.15 * (ymax - ymin);
                      ymin -= pad;
                      ymax += pad;
                      if (ymin < 0.0) ymin = 0.0;
                      if (ymin > 1.0) ymin = 0.9;
                      if (ymax < 1.0) ymax = 1.1;
                    }

                    // Save the r04 half-closure y-range to be reused by photon half-closure plots
                    if (rKey == "r04")
                    {
                      gYmin_r04_halfClosure = ymin;
                      gYmax_r04_halfClosure = ymax;
                    }

                    TCanvas c(TString::Format("c_halfClosure_vs_pt_%s", rKey.c_str()).Data(), "c_halfClosure_vs_pt", 900, 700);
                    ApplyCanvasMargins1D(c);

                    TGraphErrors g((int)xPt.size(), &xPt[0], &yRat[0], &exPt[0], &eyRat[0]);
                    g.SetLineWidth(2);
                    g.SetMarkerStyle(20);
                    g.SetTitle("");
                    g.Draw("AP");

                    g.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV] (bin centers)");
                    g.GetXaxis()->SetTitleOffset(1.15);
                    g.GetYaxis()->SetTitle("Half-closure: Unfolded MC / Truth MC");
                    g.GetYaxis()->SetRangeUser(ymin, ymax);

                    TLine l1(g.GetXaxis()->GetXmin(), 1.0, g.GetXaxis()->GetXmax(), 1.0);
                    l1.SetLineStyle(2);
                    l1.SetLineWidth(2);
                    l1.Draw("same");

                    // Centered title
                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(22);
                      tx.SetTextSize(0.040);
                      tx.DrawLatex(0.50, 0.965,
                                   TString::Format("Half-closure test, train(A) unfold(B), R = %.1f", R).Data());
                    }

                    // Top-left Bayes annotation
                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(13);
                      tx.SetTextSize(0.038);
                      tx.DrawLatex(0.20, 0.90, TString::Format("Bayes it=%d (xJ)", kBayesIterXJ).Data());
                      tx.DrawLatex(0.20, 0.855, "2D (p_{T}^{#gamma}, x_{J}) unfolding");
                    }

                    SaveCanvas(c, JoinPath(closurePlotsOut, "halfClosure_unfoldedOverTruth_integral_vs_pTgamma.png"));
                }

                if (h2Truth_halfB) delete h2Truth_halfB;
                if (h2UnfoldTruth_half) delete h2UnfoldTruth_half;
                if (hUnfoldTruthGlob_half) delete hUnfoldTruthGlob_half;
                if (hMeasSimGlob_halfB_meas) delete hMeasSimGlob_halfB_meas;
              }

              if (hRsp_measXtruth_A) delete hRsp_measXtruth_A;
              if (hRsp_measXtruth_B) delete hRsp_measXtruth_B;
              if (hMeasSimGlob_A) delete hMeasSimGlob_A;
              if (hTruthSimGlob_A) delete hTruthSimGlob_A;
              if (hMeasSimGlob_B) delete hMeasSimGlob_B;
              if (hTruthSimGlob_B) delete hTruthSimGlob_B;
            }
        }

        for (int i = 0; i < nPtAll; ++i)
        {
          const PtBin& b = analysisRecoBins[i];
          const double cen = 0.5 * (b.lo + b.hi);

          if (!h2UnfoldTruth || !hPhoUnfoldTruth) continue;

          const int ixTruth = h2UnfoldTruth->GetXaxis()->FindBin(cen);
          const int ibPho   = hPhoUnfoldTruth->GetXaxis()->FindBin(cen);

          if (ixTruth < 1 || ixTruth > h2UnfoldTruth->GetXaxis()->GetNbins()) continue;
          if (ibPho   < 1 || ibPho   > hPhoUnfoldTruth->GetXaxis()->GetNbins()) continue;

          const double Npho  = hPhoUnfoldTruth->GetBinContent(ibPho);
          const double eNpho = hPhoUnfoldTruth->GetBinError  (ibPho);

          const string labCanon = TString::Format("%d-%d GeV", b.lo, b.hi).Data();
          const string labTruth = AxisBinLabel(h2UnfoldTruth->GetXaxis(), ixTruth, "GeV", 0);
          const string labPho   = AxisBinLabel(hPhoUnfoldTruth->GetXaxis(), ibPho, "GeV", 0);

          TH1D* hXJ = h2UnfoldTruth->ProjectionY(
            TString::Format("h_xJ_unfTruth_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
            ixTruth, ixTruth, "e"
          );
          if (!hXJ) continue;
          hXJ->SetDirectory(nullptr);
          EnsureSumw2(hXJ);

          TH1D* hPerPho = (TH1D*)hXJ->Clone(
            TString::Format("h_xJ_unf_perPho_%s_pTbin%d", rKey.c_str(), i + 1).Data()
          );
          hPerPho->SetDirectory(nullptr);
          EnsureSumw2(hPerPho);

          for (int ib = 0; ib <= hPerPho->GetNbinsX() + 1; ++ib)
          {
            if (ib == 0 || ib == hPerPho->GetNbinsX() + 1)
            {
              hPerPho->SetBinContent(ib, 0.0);
              hPerPho->SetBinError  (ib, 0.0);
              continue;
            }

            const double num  = hXJ->GetBinContent(ib);
            const double eNum = hXJ->GetBinError  (ib);
            const double wid  = hXJ->GetBinWidth  (ib);

            if (Npho <= 0.0 || wid <= 0.0)
            {
              hPerPho->SetBinContent(ib, 0.0);
              hPerPho->SetBinError  (ib, 0.0);
              continue;
            }

            const double val = num / (Npho * wid);

            double relNum = 0.0;
            if (num > 0.0) relNum = eNum / num;

            double relDen = 0.0;
            if (Npho > 0.0) relDen = eNpho / Npho;

            const double err = val * std::sqrt(relNum*relNum + relDen*relDen);

            hPerPho->SetBinContent(ib, val);
            hPerPho->SetBinError  (ib, err);
          }

            hPerPho->SetTitle("");
            hPerPho->GetXaxis()->SetTitle("x_{J}");
            hPerPho->GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
            hPerPho->SetLineWidth(2);
            hPerPho->SetMarkerStyle(20);
            hPerPho->SetMarkerSize(0.85);

            perPhoHists[i] = hPerPho;

            if (gApplyPurityCorrectionForUnfolding && h2UnfoldTruth_jetEffCorr)
            {
              TH1D* hXJ_jetEffCorr = h2UnfoldTruth_jetEffCorr->ProjectionY(
                TString::Format("h_xJ_unfTruth_jetEffCorr_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
                ixTruth, ixTruth, "e"
              );

              if (hXJ_jetEffCorr)
              {
                hXJ_jetEffCorr->SetDirectory(nullptr);
                EnsureSumw2(hXJ_jetEffCorr);

                TH1D* hPerPhoJetEffCorr = (TH1D*)hXJ_jetEffCorr->Clone(
                  TString::Format("h_xJ_unf_perPho_jetEffCorr_%s_pTbin%d", rKey.c_str(), i + 1).Data()
                );

                if (hPerPhoJetEffCorr)
                {
                  hPerPhoJetEffCorr->SetDirectory(nullptr);
                  EnsureSumw2(hPerPhoJetEffCorr);

                  for (int ib = 0; ib <= hPerPhoJetEffCorr->GetNbinsX() + 1; ++ib)
                  {
                    if (ib == 0 || ib == hPerPhoJetEffCorr->GetNbinsX() + 1)
                    {
                      hPerPhoJetEffCorr->SetBinContent(ib, 0.0);
                      hPerPhoJetEffCorr->SetBinError  (ib, 0.0);
                      continue;
                    }

                    const double num  = hXJ_jetEffCorr->GetBinContent(ib);
                    const double eNum = hXJ_jetEffCorr->GetBinError  (ib);
                    const double wid  = hXJ_jetEffCorr->GetBinWidth  (ib);

                    if (!(std::isfinite(num) && std::isfinite(eNum) &&
                          std::isfinite(Npho) && std::isfinite(eNpho) &&
                          Npho > 0.0 && wid > 0.0))
                    {
                      hPerPhoJetEffCorr->SetBinContent(ib, 0.0);
                      hPerPhoJetEffCorr->SetBinError  (ib, 0.0);
                      continue;
                    }

                    const double val = num / (Npho * wid);
                    const double var = (eNum * eNum) / (Npho * Npho * wid * wid)
                                     + (num * num * eNpho * eNpho) / (Npho * Npho * Npho * Npho * wid * wid);
                    const double err = (var > 0.0 && std::isfinite(var)) ? std::sqrt(var) : 0.0;

                    hPerPhoJetEffCorr->SetBinContent(ib, val);
                    hPerPhoJetEffCorr->SetBinError  (ib, err);
                  }

                  hPerPhoJetEffCorr->SetTitle("");
                  hPerPhoJetEffCorr->GetXaxis()->SetTitle("x_{J}");
                  hPerPhoJetEffCorr->GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                  hPerPhoJetEffCorr->SetLineWidth(2);
                  hPerPhoJetEffCorr->SetMarkerStyle(24);
                  hPerPhoJetEffCorr->SetMarkerSize(0.85);
                  hPerPhoJetEffCorr->SetMarkerColor(kRed + 1);
                  hPerPhoJetEffCorr->SetLineColor(kRed + 1);

                  perPhoHists_jetEffCorr[i] = hPerPhoJetEffCorr;
                }

                delete hXJ_jetEffCorr;
              }
            }

          cout << ANSI_BOLD_YEL
                 << "[PER-PHOTON DEBUG] rKey=" << rKey
                 << "  storedPtIndex=" << i
                 << "  canonicalRecoPt=" << b.lo << "-" << b.hi
                 << "  truthXbin=" << ixTruth
                 << "  phoTruthBin=" << ibPho
                 << "  Npho=" << Npho
                 << "  eNpho=" << eNpho
                 << "  rawIntegral=" << hXJ->Integral(1, hXJ->GetNbinsX())
                 << "  rawIntegral(width)=" << hXJ->Integral(1, hXJ->GetNbinsX(), "width")
                 << "  perPhoIntegral(width)=" << hPerPho->Integral(1, hPerPho->GetNbinsX(), "width")
                 << ANSI_RESET << "\n";

          DumpTH1Summary(
              TString::Format("Raw unfolded xJ before per-photon normalization (%s, pT %d-%d)", rKey.c_str(), b.lo, b.hi).Data(),
              hXJ
            );
          DumpTH1Summary(
              TString::Format("Per-photon unfolded xJ used for plots (%s, pT %d-%d)", rKey.c_str(), b.lo, b.hi).Data(),
              hPerPho
            );

          // -------------------------------------------------------------------
          //  covariance-error version of the same unfolded per-photon spectrum
          //      (used only for ToyUnfoldingVsCovariance overlays)
          // -------------------------------------------------------------------
          if (h2UnfoldTruth_cov && hPhoUnfoldTruth_cov)
          {
              const int ixTruthCov = h2UnfoldTruth_cov->GetXaxis()->FindBin(cen);
              const int ibPhoCov   = hPhoUnfoldTruth_cov->GetXaxis()->FindBin(cen);

              if (ixTruthCov >= 1 && ixTruthCov <= h2UnfoldTruth_cov->GetXaxis()->GetNbins() &&
                  ibPhoCov   >= 1 && ibPhoCov   <= hPhoUnfoldTruth_cov->GetXaxis()->GetNbins())
              {
                const double NphoCov  = hPhoUnfoldTruth_cov->GetBinContent(ibPhoCov);
                const double eNphoCov = hPhoUnfoldTruth_cov->GetBinError  (ibPhoCov);

                TH1D* hXJ_cov = h2UnfoldTruth_cov->ProjectionY(
                  TString::Format("h_xJ_unfTruth_cov_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
                  ixTruthCov, ixTruthCov, "e"
                );
                if (hXJ_cov)
                {
                  hXJ_cov->SetDirectory(nullptr);
                  EnsureSumw2(hXJ_cov);

                  TH1D* hPerPhoCov = (TH1D*)hXJ_cov->Clone(
                    TString::Format("h_xJ_unf_perPho_cov_%s_pTbin%d", rKey.c_str(), i + 1).Data()
                  );
                  hPerPhoCov->SetDirectory(nullptr);
                  EnsureSumw2(hPerPhoCov);

                  for (int ib = 0; ib <= hPerPhoCov->GetNbinsX() + 1; ++ib)
                  {
                    if (ib == 0 || ib == hPerPhoCov->GetNbinsX() + 1)
                    {
                      hPerPhoCov->SetBinContent(ib, 0.0);
                      hPerPhoCov->SetBinError  (ib, 0.0);
                      continue;
                    }

                    const double num  = hXJ_cov->GetBinContent(ib);
                    const double eNum = hXJ_cov->GetBinError  (ib);
                    const double wid  = hXJ_cov->GetBinWidth  (ib);

                    if (NphoCov <= 0.0 || wid <= 0.0)
                    {
                      hPerPhoCov->SetBinContent(ib, 0.0);
                      hPerPhoCov->SetBinError  (ib, 0.0);
                      continue;
                    }

                    const double val = num / (NphoCov * wid);

                    double relNum = 0.0;
                    if (num > 0.0) relNum = eNum / num;

                    double relDen = 0.0;
                    if (NphoCov > 0.0) relDen = eNphoCov / NphoCov;

                    const double err = val * std::sqrt(relNum*relNum + relDen*relDen);

                    hPerPhoCov->SetBinContent(ib, val);
                    hPerPhoCov->SetBinError  (ib, err);
                  }

                  hPerPhoCov->SetTitle("");
                  hPerPhoCov->GetXaxis()->SetTitle("x_{J}");
                  hPerPhoCov->GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                  hPerPhoCov->SetLineWidth(2);
                  hPerPhoCov->SetMarkerStyle(24);
                  hPerPhoCov->SetMarkerSize(0.85);
                  hPerPhoCov->SetMarkerColor(kRed + 1);
                  hPerPhoCov->SetLineColor(kRed + 1);

                  perPhoHists_cov[i] = hPerPhoCov;

                  delete hXJ_cov;
                }
             }
          }

          // -------------------------------------------------------------------
          // "before unfolding data" (measured reco per-photon) for this pT bin
          // -------------------------------------------------------------------
          if (h2RecoData && hPhoRecoData)
          {
              const int ixReco = h2RecoData->GetXaxis()->FindBin(cen);
              const int ibPhoR = hPhoRecoData->GetXaxis()->FindBin(cen);

              if (ixReco >= 1 && ixReco <= h2RecoData->GetXaxis()->GetNbins() &&
                  ibPhoR >= 1 && ibPhoR <= hPhoRecoData->GetXaxis()->GetNbins())
              {
                const double NphoReco  = hPhoRecoData->GetBinContent(ibPhoR);
                const double eNphoReco = hPhoRecoData->GetBinError  (ibPhoR);

                TH1D* hXJ_reco = h2RecoData->ProjectionY(
                  TString::Format("h_xJ_recoData_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
                  ixReco, ixReco, "e"
                );
                if (hXJ_reco)
                {
                  hXJ_reco->SetDirectory(nullptr);
                  EnsureSumw2(hXJ_reco);

                  TH1D* hPerPhoReco = (TH1D*)hXJ_reco->Clone(
                    TString::Format("h_xJ_recoData_perPho_%s_pTbin%d", rKey.c_str(), i + 1).Data()
                  );
                  hPerPhoReco->SetDirectory(nullptr);
                  EnsureSumw2(hPerPhoReco);

                  for (int ib = 0; ib <= hPerPhoReco->GetNbinsX() + 1; ++ib)
                  {
                    if (ib == 0 || ib == hPerPhoReco->GetNbinsX() + 1)
                    {
                      hPerPhoReco->SetBinContent(ib, 0.0);
                      hPerPhoReco->SetBinError  (ib, 0.0);
                      continue;
                    }

                    const double num  = hXJ_reco->GetBinContent(ib);
                    const double eNum = hXJ_reco->GetBinError  (ib);
                    const double wid  = hXJ_reco->GetBinWidth  (ib);

                    if (NphoReco <= 0.0 || wid <= 0.0)
                    {
                      hPerPhoReco->SetBinContent(ib, 0.0);
                      hPerPhoReco->SetBinError  (ib, 0.0);
                      continue;
                    }

                    const double val = num / (NphoReco * wid);

                    double relNum = 0.0;
                    if (num > 0.0) relNum = eNum / num;

                    double relDen = 0.0;
                    if (NphoReco > 0.0) relDen = eNphoReco / NphoReco;

                    const double err = val * std::sqrt(relNum*relNum + relDen*relDen);

                    hPerPhoReco->SetBinContent(ib, val);
                    hPerPhoReco->SetBinError  (ib, err);
                  }

                  hPerPhoReco->SetTitle("");
                  hPerPhoReco->GetXaxis()->SetTitle("x_{J}");
                  hPerPhoReco->GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                  hPerPhoReco->SetLineWidth(2);
                  hPerPhoReco->SetMarkerStyle(24);
                  hPerPhoReco->SetMarkerSize(0.85);
                  hPerPhoReco->SetMarkerColor(kBlue + 1);
                  hPerPhoReco->SetLineColor(kBlue + 1);

                  perPhoBeforeDataHists[i] = hPerPhoReco;

                  delete hXJ_reco;
                }
              }
            }

            // -------------------------------------------------------------------
            // "truth MC" (SIM truth per-photon) for this pT bin
            // -------------------------------------------------------------------
            if (h2TruthSim && hPhoTruthSim)
            {
              const int ixT = h2TruthSim->GetXaxis()->FindBin(cen);
              const int ibPhoT = hPhoTruthSim->GetXaxis()->FindBin(cen);

              if (ixT >= 1 && ixT <= h2TruthSim->GetXaxis()->GetNbins() &&
                  ibPhoT >= 1 && ibPhoT <= hPhoTruthSim->GetXaxis()->GetNbins())
              {
                const double NphoT  = hPhoTruthSim->GetBinContent(ibPhoT);
                const double eNphoT = hPhoTruthSim->GetBinError  (ibPhoT);

                TH1D* hXJ_truth = h2TruthSim->ProjectionY(
                  TString::Format("h_xJ_truthSim_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
                  ixT, ixT, "e"
                );
                if (hXJ_truth)
                {
                  hXJ_truth->SetDirectory(nullptr);
                  EnsureSumw2(hXJ_truth);

                  TH1D* hPerPhoTruth = (TH1D*)hXJ_truth->Clone(
                    TString::Format("h_xJ_truthSim_perPho_%s_pTbin%d", rKey.c_str(), i + 1).Data()
                  );
                  hPerPhoTruth->SetDirectory(nullptr);
                  EnsureSumw2(hPerPhoTruth);

                  for (int ib = 0; ib <= hPerPhoTruth->GetNbinsX() + 1; ++ib)
                  {
                    if (ib == 0 || ib == hPerPhoTruth->GetNbinsX() + 1)
                    {
                      hPerPhoTruth->SetBinContent(ib, 0.0);
                      hPerPhoTruth->SetBinError  (ib, 0.0);
                      continue;
                    }

                    const double num  = hXJ_truth->GetBinContent(ib);
                    const double eNum = hXJ_truth->GetBinError  (ib);
                    const double wid  = hXJ_truth->GetBinWidth  (ib);

                    if (NphoT <= 0.0 || wid <= 0.0)
                    {
                      hPerPhoTruth->SetBinContent(ib, 0.0);
                      hPerPhoTruth->SetBinError  (ib, 0.0);
                      continue;
                    }

                    const double val = num / (NphoT * wid);

                    double relNum = 0.0;
                    if (num > 0.0) relNum = eNum / num;

                    double relDen = 0.0;
                    if (NphoT > 0.0) relDen = eNphoT / NphoT;

                    const double err = val * std::sqrt(relNum*relNum + relDen*relDen);

                    hPerPhoTruth->SetBinContent(ib, val);
                    hPerPhoTruth->SetBinError  (ib, err);
                  }

                  hPerPhoTruth->SetTitle("");
                  hPerPhoTruth->GetXaxis()->SetTitle("x_{J}");
                  hPerPhoTruth->GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                  hPerPhoTruth->SetLineWidth(2);
                  hPerPhoTruth->SetMarkerStyle(25);
                  hPerPhoTruth->SetMarkerSize(0.85);
                  hPerPhoTruth->SetMarkerColor(kMagenta + 1);
                  hPerPhoTruth->SetLineColor(kMagenta + 1);

                  perPhoTruthHists[i] = hPerPhoTruth;

                  delete hXJ_truth;
                }
              }
            }

            const double intNum = hXJ->Integral(1, hXJ->GetNbinsX(), "width");
            const double intPerPho = (Npho > 0.0 ? intNum / Npho : 0.0);

            lines.push_back(
                TString::Format(
                  "pT^gamma reco=%s  -> truthBin=%s  (phoTruthBin=%s)  Npho=%.6g±%.3g  Int[dN/dxJ]=%.6g  Int[(1/Npho)dN/dxJ]=%.6g",
                  labCanon.c_str(), labTruth.c_str(), labPho.c_str(),
                  Npho, eNpho,
                  intNum, intPerPho
                ).Data()
            );

            {
                TCanvas c(TString::Format("c_perPho_%s_%d", rKey.c_str(), i + 1).Data(), "c_perPho", 900, 700);
                ApplyCanvasMargins1D(c);

                // Match the table-pad styling: tight Y-range from content+error, same title + TLatex block
                double maxY = 0.0;
                const int nxb = hPerPho->GetNbinsX();
                for (int ib = 1; ib <= nxb; ++ib)
                {
                  const double y  = hPerPho->GetBinContent(ib);
                  const double ey = hPerPho->GetBinError(ib);
                  const double v  = y + ey;
                  if (v > maxY) maxY = v;
                }

                hPerPho->SetMinimum(0.0);
                hPerPho->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                hPerPho->GetXaxis()->SetRangeUser(0.0, 2.0);
                hPerPho->Draw("E1");

                // Centered title (same as table pads)
                {
                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(22);
                  tx.SetTextSize(0.042);

                  tx.DrawLatex(0.52, 0.955,
                               TString::Format("Per-photon particle-level x_{J#gamma}, p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                               b.lo, b.hi, R).Data());
                }

                // Per-pad cut/trigger label (top-right, 3 lines) (same as table pads)
                {
                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(31);
                  tx.SetTextSize(0.04);

                    const double xR = 0.93;
                    tx.DrawLatex(xR, 0.60, TString::Format("Bayes it = %d", kBayesIterXJ).Data());
                    tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                    tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                    tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                    tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                }

                SaveCanvas(c, JoinPath(finalUnfoldedOut, TString::Format("xJ_unfolded_perPhoton_pTbin%d.png", i + 1).Data()));

                // -------------------------------------------------------------------
                //  Toy unfolding vs analytic covariance errors (DATA)
                //   output: <rOut>/ToyUnfoldingVsCovariance/xJ_ToyUnfoldingVsCovariance_1x2_pTbin%d.png
                //
                //   Show as a 1x2 canvas:
                //     LHS = analytic covariance (kCovariance)
                //     RHS = toy unfolding errors (kCovToy)
                //
                //  Additionally output the ratio of statistical error magnitudes vs xJ:
                //     ratio(xJ) = err_cov(xJ) / err_toy(xJ)
                //   output: <rOut>/ToyUnfoldingVsCovariance/xJ_ToyUnfoldingVsCovariance_errorRatio_pTbin%d.png
                // -------------------------------------------------------------------
                if (perPhoHists[i] && perPhoHists_cov[i])
                {
                    TCanvas cTC(
                      TString::Format("c_toyVsCov_1x2_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
                      "c_toyVsCov_1x2", 1800, 700
                    );
                    cTC.Divide(2, 1, 0.001, 0.001);

                    TH1* hToy = (TH1*)perPhoHists[i]->Clone(
                      TString::Format("%s_clone_toyVsCov_1x2_toy_pTbin%d", perPhoHists[i]->GetName(), i + 1).Data()
                    );
                    TH1* hCov = (TH1*)perPhoHists_cov[i]->Clone(
                      TString::Format("%s_clone_toyVsCov_1x2_cov_pTbin%d", perPhoHists_cov[i]->GetName(), i + 1).Data()
                    );

                    if (hToy) { hToy->SetDirectory(nullptr); EnsureSumw2(hToy); }
                    if (hCov) { hCov->SetDirectory(nullptr); EnsureSumw2(hCov); }

                    if (hToy && hCov)
                    {
                      hToy->GetXaxis()->SetRangeUser(0.0, 2.0);
                      hCov->GetXaxis()->SetRangeUser(0.0, 2.0);

                      double maxY = 0.0;
                      const int nxb = hToy->GetNbinsX();
                      for (int ib = 1; ib <= nxb; ++ib)
                      {
                        const double y1  = hToy->GetBinContent(ib);
                        const double ey1 = hToy->GetBinError  (ib);
                        const double y2  = hCov->GetBinContent(ib);
                        const double ey2 = hCov->GetBinError  (ib);
                        if (y1 + ey1 > maxY) maxY = y1 + ey1;
                        if (y2 + ey2 > maxY) maxY = y2 + ey2;
                      }
                      const double yMaxUse = (maxY > 0.0) ? (1.15 * maxY) : 1.0;

                      // -----------------------------
                      // Pad 1 (LHS): kCovariance
                      // -----------------------------
                      cTC.cd(1);
                      gPad->SetLeftMargin(0.12);
                      gPad->SetRightMargin(0.04);
                      gPad->SetBottomMargin(0.12);
                      gPad->SetTopMargin(0.06);

                      hCov->SetMinimum(0.0);
                      hCov->SetMaximum(yMaxUse);
                      hCov->Draw("E1");

                      // Centered title
                      {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(22);
                        tx.SetTextSize(0.040);
                        tx.DrawLatex(0.50, 0.965,
                                     TString::Format("Per-photon particle-level x_{J#gamma}, p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                                     b.lo, b.hi, R).Data());
                      }

                      // Method label (top-left)
                      {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(13);
                        tx.SetTextSize(0.036);
                        tx.DrawLatex(0.15, 0.88, "Analytic covariance (kCovariance)");
                      }

                      // Per-pad cut/trigger label (top-right)
                      {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(31);
                        tx.SetTextSize(0.035);

                        const double xR = 0.92;
                        tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                        tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                        tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                        tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                        tx.DrawLatex(xR, 0.60, TString::Format("Bayes it = %d", kBayesIterXJ).Data());
                      }

                      // -----------------------------
                      // Pad 2 (RHS): kCovToy
                      // -----------------------------
                      cTC.cd(2);
                      gPad->SetLeftMargin(0.12);
                      gPad->SetRightMargin(0.04);
                      gPad->SetBottomMargin(0.12);
                      gPad->SetTopMargin(0.06);

                      hToy->SetMinimum(0.0);
                      hToy->SetMaximum(yMaxUse);
                      hToy->Draw("E1");

                      // Centered title
                      {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(22);
                        tx.SetTextSize(0.040);
                        tx.DrawLatex(0.50, 0.965,
                                     TString::Format("Per-photon particle-level x_{J#gamma}, p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                                     b.lo, b.hi, R).Data());
                      }

                      // Method label (top-left)
                      {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(13);
                        tx.SetTextSize(0.036);
                        tx.DrawLatex(0.15, 0.88, "Toy unfolding errors (kCovToy)");
                      }

                      // Per-pad cut/trigger label (top-right)
                      {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(31);
                        tx.SetTextSize(0.035);

                        const double xR = 0.92;
                        tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                        tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                        tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                        tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                        tx.DrawLatex(xR, 0.60, TString::Format("Bayes it = %d", kBayesIterXJ).Data());
                      }

                      SaveCanvas(cTC, JoinPath(toyVsCovOut, TString::Format("xJ_ToyUnfoldingVsCovariance_1x2_pTbin%d.png", i + 1).Data()));

                      // -----------------------------
                      // Error-magnitude ratio vs xJ: err_cov / err_toy
                      // -----------------------------
                      TH1* hRatio = (TH1*)hCov->Clone(
                        TString::Format("h_errRatio_covOverToy_%s_pTbin%d", rKey.c_str(), i + 1).Data()
                      );
                      if (hRatio)
                      {
                        hRatio->SetDirectory(nullptr);
                        EnsureSumw2(hRatio);

                          const double relToySigma =
                            (kNToysXJFinal > 1) ? (1.0 / std::sqrt(2.0 * (double)(kNToysXJFinal - 1))) : 0.0;

                          for (int ib = 0; ib <= hRatio->GetNbinsX() + 1; ++ib)
                          {
                            if (ib == 0 || ib == hRatio->GetNbinsX() + 1)
                            {
                              hRatio->SetBinContent(ib, 0.0);
                              hRatio->SetBinError  (ib, 0.0);
                              continue;
                            }

                            const double eCov = hCov->GetBinError(ib);
                            const double eToy = hToy->GetBinError(ib);

                            if (eToy > 0.0 && eCov >= 0.0)
                            {
                              const double r  = eCov / eToy;
                              const double er = std::fabs(r) * relToySigma;

                              hRatio->SetBinContent(ib, r);
                              hRatio->SetBinError  (ib, er);
                            }
                            else
                            {
                              hRatio->SetBinContent(ib, 0.0);
                              hRatio->SetBinError  (ib, 0.0);
                            }
                          }

                          hRatio->SetTitle("");
                          hRatio->GetXaxis()->SetTitle("x_{J}");
                          hRatio->GetYaxis()->SetTitle("#sigma_{cov} / #sigma_{toy}");
                          hRatio->SetLineWidth(0);
                          hRatio->SetMarkerStyle(20);
                          hRatio->SetMarkerSize(0.95);

                        // Store for optional summary table (first 6 pT bins)
                        perPhoErrRatio[i] = (TH1*)hRatio->Clone(
                          TString::Format("h_errRatio_covOverToy_store_%s_pTbin%d", rKey.c_str(), i + 1).Data()
                        );
                        if (perPhoErrRatio[i])
                        {
                          perPhoErrRatio[i]->SetDirectory(nullptr);
                          EnsureSumw2(perPhoErrRatio[i]);
                        }

                        TCanvas cR(
                          TString::Format("c_errRatio_%s_pTbin%d", rKey.c_str(), i + 1).Data(),
                          "c_errRatio", 900, 700
                        );
                        ApplyCanvasMargins1D(cR);

                        hRatio->GetXaxis()->SetRangeUser(0.0, 2.0);
                        hRatio->SetMinimum(0.0);

                        double rMax = 0.0;
                        for (int ib = 1; ib <= hRatio->GetNbinsX(); ++ib)
                        {
                          const double v = hRatio->GetBinContent(ib);
                          if (v > rMax) rMax = v;
                        }
                        hRatio->SetMaximum((rMax > 0.0) ? (1.25 * rMax) : 2.0);

                        hRatio->Draw("P E1");

                        // Centered title
                        {
                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextAlign(22);
                          tx.SetTextSize(0.040);
                          tx.DrawLatex(0.50, 0.965,
                                       TString::Format("Error ratio, p_{T}^{#gamma} %d-%d GeV, R = %.1f", b.lo, b.hi, R).Data());
                        }

                        // Per-pad cut/trigger label (top-right)
                        {
                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextAlign(31);
                          tx.SetTextSize(0.035);

                          const double xR = 0.92;
                          tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                          tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                          tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                          tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                          tx.DrawLatex(xR, 0.60, TString::Format("Bayes it = %d", kBayesIterXJ).Data());
                        }

                        SaveCanvas(cR, JoinPath(toyVsCovOut, TString::Format("xJ_ToyUnfoldingVsCovariance_errorRatio_pTbin%d.png", i + 1).Data()));

                        delete hRatio;
                      }
                    }

                    if (hToy) delete hToy;
                    if (hCov) delete hCov;
                }
                // -------------------------------------------------------------------
                // before/after unfolding overlay (DATA)
                // -------------------------------------------------------------------
                if (perPhoBeforeDataHists[i] && perPhoHists[i])
                {
                  TCanvas cBA(TString::Format("c_beforeAfterUnf_data_%s_%d", rKey.c_str(), i + 1).Data(),
                             "c_beforeAfterUnf_data", 900, 700);
                  ApplyCanvasMargins1D(cBA);

                  TH1* hA = (TH1*)perPhoBeforeDataHists[i]->Clone(
                    TString::Format("hTmp_beforeUnf_data_%s_%d", rKey.c_str(), i + 1).Data()
                  );
                  TH1* hB = (TH1*)perPhoHists[i]->Clone(
                    TString::Format("hTmp_unfolded_data_%s_%d", rKey.c_str(), i + 1).Data()
                  );

                  if (hA && hB)
                  {
                    hA->SetDirectory(nullptr);
                    hB->SetDirectory(nullptr);

                    double maxY = 0.0;
                    for (int ib = 1; ib <= hA->GetNbinsX(); ++ib)
                    {
                      const double v1 = hA->GetBinContent(ib) + hA->GetBinError(ib);
                      const double v2 = hB->GetBinContent(ib) + hB->GetBinError(ib);
                      if (v1 > maxY) maxY = v1;
                      if (v2 > maxY) maxY = v2;
                    }

                    hA->SetTitle("");
                    hA->SetMinimum(0.0);
                    hA->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                    hA->GetXaxis()->SetRangeUser(0.0, 2.0);
                    hA->Draw("E1");
                    hB->Draw("E1 same");

                    TLegend leg(0.46, 0.32, 0.92, 0.48);
                    leg.SetBorderSize(0);
                    leg.SetFillStyle(0);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.040);
                    leg.AddEntry(hA, "before unfolding data", "pe");
                    leg.AddEntry(hB, "unfolded data",        "pe");
                    leg.Draw();

                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(22);
                      tx.SetTextSize(0.042);

                      tx.DrawLatex(0.52, 0.955,
                                   TString::Format("Before vs after unfolding (DATA), p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                                   b.lo, b.hi, R).Data());
                    }

                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(31);
                      tx.SetTextSize(0.04);

                      const double xR = 0.93;
                      tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                      tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                      tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                      tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                    }

                    SaveCanvas(cBA, JoinPath(beforeAfterDataOut,
                      TString::Format("xJ_before_after_unfoldingOverlay_data_pTbin%d.png", i + 1).Data()
                    ));

                    TH1* hRatioBA = BuildRatioHist(
                      hA, hB,
                      TString::Format("h_ratio_beforeOverUnfolded_%s_pTbin%d", rKey.c_str(), i + 1).Data()
                    );

                    if (hRatioBA)
                    {
                      hRatioBA->SetDirectory(nullptr);
                      EnsureSumw2(hRatioBA);
                      hRatioBA->SetMarkerStyle(20);
                      hRatioBA->SetMarkerSize(0.95);
                      hRatioBA->SetLineWidth(2);
                      hRatioBA->SetTitle("");
                      hRatioBA->GetYaxis()->SetTitle("Before unfolding / unfolded");

                        double minR =  1e99;
                        double maxR = -1e99;
                        bool haveRatioPoint = false;

                        for (int ib = 1; ib <= hRatioBA->GetNbinsX(); ++ib)
                        {
                          const double y  = hRatioBA->GetBinContent(ib);
                          const double ey = hRatioBA->GetBinError(ib);

                          if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                          if (y <= 0.0 && ey <= 0.0) continue;

                          const double lo = y - ey;
                          const double hi = y + ey;

                          haveRatioPoint = true;
                          if (lo < minR) minR = lo;
                          if (hi > maxR) maxR = hi;
                        }

                        double yMinUse = 0.0;
                        double yMaxUse = 1.2;

                        if (haveRatioPoint)
                        {
                          yMaxUse = std::max(1.0, maxR);

                          if (!std::isfinite(yMaxUse) || yMaxUse <= 0.0)
                          {
                            yMinUse = 0.0;
                            yMaxUse = 1.2;
                          }
                          else
                          {
                            const double pad = std::max(0.08, 0.12 * yMaxUse);
                            yMinUse = 0.0;
                            yMaxUse += pad;
                          }

                          if (yMaxUse < 0.2)
                          {
                            yMaxUse = 0.2;
                          }
                        }

                      TCanvas cBAR(TString::Format("c_beforeAfterRatio_data_%s_%d", rKey.c_str(), i + 1).Data(),
                                   "c_beforeAfterRatio_data", 900, 700);
                      ApplyCanvasMargins1D(cBAR);

                      hRatioBA->SetMinimum(yMinUse);
                      hRatioBA->SetMaximum(yMaxUse);
                      hRatioBA->GetXaxis()->SetRangeUser(0.0, 2.0);
                      hRatioBA->Draw("E1");

                      TLine l1(0.0, 1.0, 2.0, 1.0);
                      l1.SetLineStyle(2);
                      l1.SetLineWidth(2);
                      l1.Draw("same");

                      {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(22);
                        tx.SetTextSize(0.042);

                        tx.DrawLatex(0.52, 0.955,
                                     TString::Format("Before / unfolded ratio, p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                                     b.lo, b.hi, R).Data());
                      }

                      {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(31);
                        tx.SetTextSize(0.04);

                        const double xR = 0.93;
                        tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                        tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                        tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                        tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                      }

                      SaveCanvas(cBAR, JoinPath(beforeAfterDataOut,
                        TString::Format("xJ_ratio_before_over_unfolded_data_pTbin%d.png", i + 1).Data()
                      ));

                      ratioBeforeVsAfterHists[i] = hRatioBA;
                    }

                    delete hA;
                    delete hB;
                  }
                }

                // -------------------------------------------------------------------
                // truth MC vs unfolded data overlay
                // -------------------------------------------------------------------
                if (perPhoTruthHists[i] && perPhoHists[i])
                {
                  TCanvas cTU(TString::Format("c_truthVsUnf_%s_%d", rKey.c_str(), i + 1).Data(),
                             "c_truthVsUnf", 900, 700);
                  ApplyCanvasMargins1D(cTU);

                  TH1* hT = (TH1*)perPhoTruthHists[i]->Clone(
                    TString::Format("hTmp_truth_%s_%d", rKey.c_str(), i + 1).Data()
                  );
                  TH1* hU = (TH1*)perPhoHists[i]->Clone(
                    TString::Format("hTmp_unfolded_data2_%s_%d", rKey.c_str(), i + 1).Data()
                  );

                  if (hT && hU)
                  {
                    hT->SetDirectory(nullptr);
                    hU->SetDirectory(nullptr);

                    double maxY = 0.0;
                    for (int ib = 1; ib <= hT->GetNbinsX(); ++ib)
                    {
                      const double v1 = hT->GetBinContent(ib) + hT->GetBinError(ib);
                      const double v2 = hU->GetBinContent(ib) + hU->GetBinError(ib);
                      if (v1 > maxY) maxY = v1;
                      if (v2 > maxY) maxY = v2;
                    }

                    hT->SetTitle("");
                    hT->SetMinimum(0.0);
                    hT->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                    hT->GetXaxis()->SetRangeUser(0.0, 2.0);
                    hT->Draw("E1");
                    hU->SetMarkerStyle(20);
                    hU->SetMarkerColor(kBlue);
                    hU->SetLineColor(kBlue);
                    hU->Draw("E1 same");

                    TLegend leg(0.69, 0.33, 0.91, 0.58);
                    leg.SetBorderSize(0);
                    leg.SetFillStyle(0);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.038);
                    leg.AddEntry(hT, "truth MC",      "pe");
                    leg.AddEntry(hU, "unfolded data", "pe");
                    leg.Draw();

                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(22);
                      tx.SetTextSize(0.042);

                      tx.DrawLatex(0.52, 0.955,
                                   TString::Format("Truth MC vs unfolded data, p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                                   b.lo, b.hi, R).Data());
                    }

                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(31);
                      tx.SetTextSize(0.04);

                      const double xR = 0.93;
                      tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                      tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                      tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                      tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                    }

                    SaveCanvas(cTU, JoinPath(beforeAfterTruthOut,
                      TString::Format("xJ_before_after_unfoldingOverlay_truth_pTbin%d.png", i + 1).Data()
                    ));

                      TH1* hRatioTU = BuildRatioHist(
                        hT, hU,
                        TString::Format("h_ratio_truthOverUnfolded_%s_pTbin%d", rKey.c_str(), i + 1).Data()
                      );

                      if (hRatioTU)
                      {
                        hRatioTU->SetDirectory(nullptr);
                        EnsureSumw2(hRatioTU);
                        hRatioTU->SetMarkerStyle(20);
                        hRatioTU->SetMarkerSize(0.95);
                        hRatioTU->SetLineWidth(2);
                        hRatioTU->SetTitle("");
                        hRatioTU->GetYaxis()->SetTitle("Truth MC / unfolded data");

                        cout << ANSI_BOLD_CYN
                             << "\n[TRUTH/UNFOLDED RATIO DEBUG] Building xJ_ratio_truthMC_over_unfolded_data_pTbin" << (i + 1)
                             << ".png"
                             << "\n  rKey=" << rKey
                             << "  R=" << R
                             << "  pTgamma=" << b.lo << "-" << b.hi << " GeV"
                             << "\n  outputDir=" << beforeAfterTruthOut
                             << "\n"
                             << ANSI_RESET;

                        cout << "  -------------------------------------------------------------------------------------------------------------\n";
                        cout << "   ib          xJ range            truth MC         errT         unfolded         errU         ratio       errR   axis\n";
                        cout << "  -------------------------------------------------------------------------------------------------------------\n";

                        double minR =  1e99;
                        double maxR = -1e99;
                        bool haveRatioPoint = false;

                        double minRRobust =  1e99;
                        double maxRRobust = -1e99;
                        bool haveRobustPoint = false;

                        for (int ib = 1; ib <= hRatioTU->GetNbinsX(); ++ib)
                        {
                            const double xLo = hRatioTU->GetXaxis()->GetBinLowEdge(ib);
                            const double xHi = hRatioTU->GetXaxis()->GetBinUpEdge(ib);

                            const double t   = hT->GetBinContent(ib);
                            const double eT  = hT->GetBinError  (ib);
                            const double u   = hU->GetBinContent(ib);
                            const double eU  = hU->GetBinError  (ib);

                            const double y   = hRatioTU->GetBinContent(ib);
                            const double ey  = hRatioTU->GetBinError  (ib);

                            bool finiteAll = std::isfinite(t) && std::isfinite(eT) &&
                                             std::isfinite(u) && std::isfinite(eU) &&
                                             std::isfinite(y) && std::isfinite(ey);

                            string axisTag = "SKIP";

                            if (finiteAll && y > 0.0 && ey > 0.0)
                            {
                              const double lo = y - ey;
                              const double hi = y + ey;

                              haveRatioPoint = true;
                              if (lo < minR) minR = lo;
                              if (hi > maxR) maxR = hi;

                              const double relErrU = (u > 0.0) ? (eU / u) : 1e99;
                              const double relErrR = (y > 0.0) ? (ey / y) : 1e99;

                              if (u > 0.0 && relErrU < 0.60 && relErrR < 0.60 && y < 5.0 && hi < 5.0)
                              {
                                axisTag = "USE";
                                haveRobustPoint = true;

                                const double loAxis = std::max(0.05, lo);
                                const double hiAxis = hi;

                                if (loAxis < minRRobust) minRRobust = loAxis;
                                if (hiAxis > maxRRobust) maxRRobust = hiAxis;
                              }
                              else if (u <= 0.0)
                              {
                                axisTag = "DEN~0";
                              }
                              else if (relErrU >= 0.60)
                              {
                                axisTag = "UNSTABLE_DEN";
                              }
                              else if (relErrR >= 0.60)
                              {
                                axisTag = "HUGE_RATIO_ERR";
                              }
                              else if (y >= 5.0 || hi >= 5.0)
                              {
                                axisTag = "OUTLIER";
                              }
                            }

                            cout << TString::Format(
                                      "  %4d   [%6.3f,%6.3f]   %12.6g %11.6g %12.6g %11.6g %12.6g %10.6g   %s\n",
                                      ib, xLo, xHi, t, eT, u, eU, y, ey, axisTag.c_str()
                                    ).Data();
                          }

                          cout << "  -------------------------------------------------------------------------------------------------------------\n";
                          cout << "  Full-axis scan:    haveRatioPoint=" << (haveRatioPoint ? "true" : "false")
                               << "  minR=" << minR
                               << "  maxR=" << maxR << "\n";
                          cout << "  Robust-axis scan:  haveRobustPoint=" << (haveRobustPoint ? "true" : "false")
                               << "  minRRobust=" << minRRobust
                               << "  maxRRobust=" << maxRRobust << "\n";

                          double yMinUse = 0.0;
                          double yMaxUse = 1.2;

                          if (haveRobustPoint)
                          {
                            yMaxUse = std::max(1.0, maxRRobust);

                            if (!std::isfinite(yMaxUse) || yMaxUse <= 0.0)
                            {
                              yMinUse = 0.0;
                              yMaxUse = 1.2;
                            }
                            else
                            {
                              const double pad = std::max(0.08, 0.12 * yMaxUse);
                              yMinUse = 0.0;
                              yMaxUse += pad;
                            }

                            if (yMaxUse < 0.20)
                            {
                              yMaxUse = 0.20;
                            }
                          }
                          else if (haveRatioPoint)
                          {
                            yMaxUse = std::max(1.0, maxR);

                            if (!std::isfinite(yMaxUse) || yMaxUse <= 0.0)
                            {
                              yMinUse = 0.0;
                              yMaxUse = 1.2;
                            }
                            else
                            {
                              const double pad = std::max(0.08, 0.12 * yMaxUse);
                              yMinUse = 0.0;
                              yMaxUse += pad;
                            }

                            if (yMaxUse < 0.20)
                            {
                              yMaxUse = 0.20;
                            }
                          }

                        cout << ANSI_BOLD_YEL
                             << "[TRUTH/UNFOLDED RATIO DEBUG] Final plotted y-range for pTbin" << (i + 1)
                             << ": yMinUse=" << yMinUse
                             << "  yMaxUse=" << yMaxUse
                             << ANSI_RESET << "\n";

                        TCanvas cTUR(TString::Format("c_truthVsUnfRatio_%s_%d", rKey.c_str(), i + 1).Data(),
                                     "c_truthVsUnfRatio", 900, 700);
                        ApplyCanvasMargins1D(cTUR);

                        hRatioTU->SetMinimum(yMinUse);
                        hRatioTU->SetMaximum(yMaxUse);
                        hRatioTU->GetXaxis()->SetRangeUser(0.0, 2.0);
                        hRatioTU->Draw("E1");

                        TLine l1(0.0, 1.0, 2.0, 1.0);
                        l1.SetLineStyle(2);
                        l1.SetLineWidth(2);
                        l1.Draw("same");

                        {
                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextAlign(22);
                          tx.SetTextSize(0.042);

                          tx.DrawLatex(0.52, 0.955,
                                       TString::Format("Truth MC / unfolded data, p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                                       b.lo, b.hi, R).Data());
                        }

                        {
                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextAlign(31);
                          tx.SetTextSize(0.04);

                          const double xR = 0.93;
                          tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                          tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                          tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                          tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                        }

                        SaveCanvas(cTUR, JoinPath(beforeAfterTruthOut,
                          TString::Format("xJ_ratio_truthMC_over_unfolded_data_pTbin%d.png", i + 1).Data()
                        ));

                        ratioTruthVsUnfoldedHists[i] = hRatioTU;
                      }

                    delete hT;
                    delete hU;
                  }
                }

                // -------------------------------------------------------------------
                // LHC overlay (ATLAS pp, HEPData ins1694678 Table 1):
                //   overlay the same ATLAS pp curve on every sPHENIX pT bin (simple first comparison)
                //   output: <rOut>/LHC_overlay/xJ_unfolded_perPhoton_LHCoverlay_pTbinX.png
                //
                //   Also make an identical jet-efficiency-corrected version when available:
                //   output: <rOut>/LHC_overlay/effCorrected/xJ_unfolded_perPhoton_LHCoverlay_pTbinX.png
                // -------------------------------------------------------------------
                if (gAtlasPP)
                {
                    const double yHeadroomScale = 1.30;
                    const double sphJetPtMin = static_cast<double>(kJetPtMin);
                    const string bbLabel = B2BLabel();
                    const double atlasJetPtMin = 31.6;

                    TCanvas cO(TString::Format("c_perPho_LHC_%s_%d", rKey.c_str(), i + 1).Data(), "c_perPho_LHC", 900, 700);
                    ApplyCanvasMargins1D(cO);
                    cO.SetLeftMargin(0.18);
                    cO.SetRightMargin(0.05);
                    cO.SetBottomMargin(0.16);
                    cO.SetTopMargin(0.08);

                    TH1* hTmp = (TH1*)hPerPho->Clone(TString::Format("hTmp_perPho_%s_%d", rKey.c_str(), i + 1).Data());
                    if (hTmp)
                    {
                      hTmp->SetDirectory(nullptr);

                      double maxY = 0.0;
                      for (int ib = 1; ib <= hTmp->GetNbinsX(); ++ib)
                      {
                        const double y  = hTmp->GetBinContent(ib);
                        const double ey = hTmp->GetBinError(ib);
                        const double v  = y + ey;
                        if (v > maxY) maxY = v;
                      }

                      for (int ip = 0; ip < gAtlasPP->GetN(); ++ip)
                      {
                        double x = 0.0, y = 0.0;
                        gAtlasPP->GetPoint(ip, x, y);
                        const double ey = gAtlasPP->GetErrorYhigh(ip);
                        const double v  = y + ey;
                        if (v > maxY) maxY = v;
                      }

                      hTmp->SetMinimum(0.0);
                      hTmp->SetMaximum((maxY > 0.0) ? (yHeadroomScale * maxY) : 1.0);
                      hTmp->GetXaxis()->SetRangeUser(0.0, 2.0);
                      hTmp->GetXaxis()->SetTitle("");
                      hTmp->GetXaxis()->SetTitleOffset(0.95);
                      hTmp->GetYaxis()->SetTitleOffset(1.10);

                      // Draw sPHENIX with horizontal error bars on the plot (TH1::Draw("E1") keeps x-errors)
                      hTmp->Draw("E1");
                      gAtlasPP->Draw("PZ same");

                        // Legend icon for sPHENIX: vertical error bar only (EX=0, EY!=0)
                        double lx = 0.0, ly = 0.0;
                        {
                          int ib0 = -1;
                          for (int ib = 1; ib <= hTmp->GetNbinsX(); ++ib)
                          {
                            if (hTmp->GetBinContent(ib) != 0.0 || hTmp->GetBinError(ib) != 0.0) { ib0 = ib; break; }
                          }
                          if (ib0 < 0) ib0 = 1;

                          lx = hTmp->GetXaxis()->GetBinCenter(ib0);
                          ly = hTmp->GetBinContent(ib0);
                        }
                        const double ley = hTmp->GetBinError(hTmp->GetXaxis()->FindBin(lx));

                        TGraphErrors gSphLeg(1);
                        gSphLeg.SetPoint(0, lx, ly);
                        gSphLeg.SetPointError(0, 0.0, ley);
                        gSphLeg.SetMarkerStyle(hTmp->GetMarkerStyle());
                        gSphLeg.SetMarkerSize(hTmp->GetMarkerSize());
                        gSphLeg.SetMarkerColor(hTmp->GetMarkerColor());
                        gSphLeg.SetLineColor(hTmp->GetLineColor());
                        gSphLeg.SetLineWidth(hTmp->GetLineWidth());

                        TLegend extraLegend(0.15, 0.78, 0.43, 0.88);
                        extraLegend.SetBorderSize(0);
                        extraLegend.SetFillStyle(0);
                        extraLegend.SetTextFont(42);
                        extraLegend.SetTextSize(0.040);
                        extraLegend.AddEntry((TObject*)nullptr, "#it{#bf{sPHENIX}} Internal", "");
                        extraLegend.AddEntry((TObject*)nullptr, "p+p #sqrt{s} = 200 GeV", "");
                        extraLegend.Draw();

                        TLegend leg(0.49,0.79,0.90,0.89);
                        leg.SetBorderSize(0);
                        leg.SetFillStyle(0);
                        leg.SetTextFont(42);
                        leg.SetTextSize(0.029);
                        leg.AddEntry(&gSphLeg,
                                       TString::Format("sPHENIX run24pp").Data(),
                                       "pe");
                        leg.AddEntry(gAtlasPP,
                                       TString::Format("ATLAS Phys. Lett. B 789 (2019) 167").Data(),
                                       "pe");
                        leg.Draw();

                        TLatex txHdr;
                        txHdr.SetNDC();
                        txHdr.SetTextFont(42);
                        txHdr.SetTextAlign(13);
                        txHdr.SetTextSize(0.040);
                        txHdr.SetTextColor(kBlack);

                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(13);
                        tx.SetTextSize(0.024);
                        tx.SetTextColor(kBlack);

                        tx.DrawLatex(0.15, 0.70, TString::Format("p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                        tx.DrawLatex(0.15, 0.65, TString::Format("R = %.1f", R).Data());
                        tx.DrawLatex(0.15, 0.60, TString::Format("|#Delta#phi| > %s", bbLabel.c_str()).Data());
                        tx.DrawLatex(0.15, 0.55, TString::Format("p_{T}^{jet} > %.0f GeV", sphJetPtMin).Data());
                        tx.DrawLatex(0.15, 0.50, TString::Format("|v_{z}| < %.0f cm", std::fabs(vzCutCm)).Data());

                        txHdr.DrawLatex(0.66, 0.72, "#it{#bf{ATLAS}}");
                        tx.DrawLatex(0.66, 0.67, "pp #sqrt{s} = 5.02 TeV");
                        tx.DrawLatex(0.66, 0.62, TString::Format("p_{T}^{#gamma} = %s GeV", kAtlasTable1PhoPtLabel.c_str()).Data());
                        tx.DrawLatex(0.66, 0.57, "R = 0.4");
                        tx.DrawLatex(0.66, 0.52, "|#Delta#phi| > 7#pi/8");
                        tx.DrawLatex(0.66, 0.47, TString::Format("p_{T}^{jet} > %.1f GeV", atlasJetPtMin).Data());

                        TLatex txAxis;
                        txAxis.SetNDC();
                        txAxis.SetTextFont(42);
                        txAxis.SetTextAlign(31);
                        txAxis.SetTextSize(0.050);
                        txAxis.DrawLatex(0.98, 0.03, "x_{J}");

                        SaveCanvas(cO, JoinPath(overlayOut, TString::Format("xJ_unfolded_perPhoton_LHCoverlay_pTbin%d.png", i + 1).Data()));

                        delete hTmp;
                  }

                  if (gApplyPurityCorrectionForUnfolding &&
                      i >= 0 &&
                      i < (int)perPhoHists_jetEffCorr.size() &&
                      perPhoHists_jetEffCorr[i])
                  {
                    const string overlayOutEffCorrected = JoinPath(overlayOut, "effCorrected");
                    EnsureDir(overlayOutEffCorrected);

                    TCanvas cOEff(TString::Format("c_perPho_LHC_effCorrected_%s_%d", rKey.c_str(), i + 1).Data(), "c_perPho_LHC_effCorrected", 900, 700);
                    ApplyCanvasMargins1D(cOEff);
                    cOEff.SetLeftMargin(0.18);
                    cOEff.SetRightMargin(0.05);
                    cOEff.SetBottomMargin(0.16);
                    cOEff.SetTopMargin(0.08);

                    TH1* hTmpEff = (TH1*)perPhoHists_jetEffCorr[i]->Clone(TString::Format("hTmp_perPho_effCorrected_%s_%d", rKey.c_str(), i + 1).Data());
                    if (hTmpEff)
                    {
                      hTmpEff->SetDirectory(nullptr);
                      hTmpEff->SetMarkerStyle(hPerPho->GetMarkerStyle());
                      hTmpEff->SetMarkerSize(hPerPho->GetMarkerSize());
                      hTmpEff->SetMarkerColor(hPerPho->GetMarkerColor());
                      hTmpEff->SetLineColor(hPerPho->GetLineColor());
                      hTmpEff->SetLineWidth(hPerPho->GetLineWidth());

                      double maxY = 0.0;
                      for (int ib = 1; ib <= hTmpEff->GetNbinsX(); ++ib)
                      {
                        const double y  = hTmpEff->GetBinContent(ib);
                        const double ey = hTmpEff->GetBinError(ib);
                        const double v  = y + ey;
                        if (v > maxY) maxY = v;
                      }

                      for (int ip = 0; ip < gAtlasPP->GetN(); ++ip)
                      {
                        double x = 0.0, y = 0.0;
                        gAtlasPP->GetPoint(ip, x, y);
                        const double ey = gAtlasPP->GetErrorYhigh(ip);
                        const double v  = y + ey;
                        if (v > maxY) maxY = v;
                      }

                      hTmpEff->SetMinimum(0.0);
                      hTmpEff->SetMaximum((maxY > 0.0) ? (yHeadroomScale * maxY) : 1.0);
                      hTmpEff->GetXaxis()->SetRangeUser(0.0, 2.0);
                      hTmpEff->GetXaxis()->SetTitle("");
                      hTmpEff->GetXaxis()->SetTitleOffset(0.95);
                      hTmpEff->GetYaxis()->SetTitleOffset(1.10);

                      // Draw sPHENIX with horizontal error bars on the plot (TH1::Draw("E1") keeps x-errors)
                      hTmpEff->Draw("E1");
                      gAtlasPP->Draw("PZ same");

                        // Legend icon for sPHENIX: vertical error bar only (EX=0, EY!=0)
                        double lx = 0.0, ly = 0.0;
                        {
                          int ib0 = -1;
                          for (int ib = 1; ib <= hTmpEff->GetNbinsX(); ++ib)
                          {
                            if (hTmpEff->GetBinContent(ib) != 0.0 || hTmpEff->GetBinError(ib) != 0.0) { ib0 = ib; break; }
                          }
                          if (ib0 < 0) ib0 = 1;

                          lx = hTmpEff->GetXaxis()->GetBinCenter(ib0);
                          ly = hTmpEff->GetBinContent(ib0);
                        }
                        const double ley = hTmpEff->GetBinError(hTmpEff->GetXaxis()->FindBin(lx));

                        TGraphErrors gSphLeg(1);
                        gSphLeg.SetPoint(0, lx, ly);
                        gSphLeg.SetPointError(0, 0.0, ley);
                        gSphLeg.SetMarkerStyle(hTmpEff->GetMarkerStyle());
                        gSphLeg.SetMarkerSize(hTmpEff->GetMarkerSize());
                        gSphLeg.SetMarkerColor(hTmpEff->GetMarkerColor());
                        gSphLeg.SetLineColor(hTmpEff->GetLineColor());
                        gSphLeg.SetLineWidth(hTmpEff->GetLineWidth());

                        TLegend extraLegend(0.15, 0.78, 0.43, 0.88);
                        extraLegend.SetBorderSize(0);
                        extraLegend.SetFillStyle(0);
                        extraLegend.SetTextFont(42);
                        extraLegend.SetTextSize(0.040);
                        extraLegend.AddEntry((TObject*)nullptr, "#it{#bf{sPHENIX}} Internal", "");
                        extraLegend.AddEntry((TObject*)nullptr, "p+p #sqrt{s} = 200 GeV", "");
                        extraLegend.Draw();

                        TLegend leg(0.49,0.79,0.90,0.89);
                        leg.SetBorderSize(0);
                        leg.SetFillStyle(0);
                        leg.SetTextFont(42);
                        leg.SetTextSize(0.029);
                        leg.AddEntry(&gSphLeg,
                                       TString::Format("sPHENIX run24pp").Data(),
                                       "pe");
                        leg.AddEntry(gAtlasPP,
                                       TString::Format("ATLAS Phys. Lett. B 789 (2019) 167").Data(),
                                       "pe");
                        leg.Draw();

                        const double sphR = RFromKey(rKey);

                        TLatex txHdr;
                        txHdr.SetNDC();
                        txHdr.SetTextFont(42);
                        txHdr.SetTextAlign(13);
                        txHdr.SetTextSize(0.040);
                        txHdr.SetTextColor(kBlack);

                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(13);
                        tx.SetTextSize(0.024);
                        tx.SetTextColor(kBlack);

                        tx.DrawLatex(0.15, 0.70, TString::Format("p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                        tx.DrawLatex(0.15, 0.65, TString::Format("R = %.1f", sphR).Data());
                        tx.DrawLatex(0.15, 0.60, TString::Format("|#Delta#phi| > %s", bbLabel.c_str()).Data());
                        tx.DrawLatex(0.15, 0.55, TString::Format("p_{T}^{jet} > %.0f GeV", sphJetPtMin).Data());
                        tx.DrawLatex(0.15, 0.50, TString::Format("|v_{z}| < %.0f cm", std::fabs(vzCutCm)).Data());

                        txHdr.DrawLatex(0.66, 0.72, "#it{#bf{ATLAS}}");
                        tx.DrawLatex(0.66, 0.67, "pp #sqrt{s} = 5.02 TeV");
                        tx.DrawLatex(0.66, 0.62, TString::Format("p_{T}^{#gamma} = %s GeV", kAtlasTable1PhoPtLabel.c_str()).Data());
                        tx.DrawLatex(0.66, 0.57, "R = 0.4");
                        tx.DrawLatex(0.66, 0.52, "|#Delta#phi| > 7#pi/8");
                        tx.DrawLatex(0.66, 0.47, TString::Format("p_{T}^{jet} > %.1f GeV", atlasJetPtMin).Data());

                        TLatex txAxis;
                        txAxis.SetNDC();
                        txAxis.SetTextFont(42);
                        txAxis.SetTextAlign(31);
                        txAxis.SetTextSize(0.050);
                        txAxis.DrawLatex(0.98, 0.03, "x_{J}");

                        SaveCanvas(cOEff, JoinPath(overlayOutEffCorrected, TString::Format("xJ_unfolded_perPhoton_LHCoverlay_pTbin%d.png", i + 1).Data()));

                        delete hTmpEff;
                    }
                  }
              }
          }

        delete hXJ;
        }

        {
              // -------------------------------------------------------------------
              // 2x3 table (first 6 pT bins): ratio of statistical error magnitudes
              //   ratio(xJ) = err_cov(xJ) / err_toy(xJ)
              //   output: <rOut>/ToyUnfoldingVsCovariance/table2x3_ToyUnfoldingVsCovariance_errorRatio.png
              // -------------------------------------------------------------------
              {
                bool anyR = false;
                const int nPadsR = nPtAll;
                for (int ii = 0; ii < nPadsR; ++ii)
                {
                    if (ii >= 0 && ii < nPtAll && perPhoErrRatio[ii]) { anyR = true; break; }
                }

                if (anyR)
                {
                  TCanvas cER(
                      TString::Format("c_tbl_errRatio_%s", rKey.c_str()).Data(),
                      "c_tbl_errRatio", 1800, 1300
                    );
                  cER.Divide(3, 3, 0.001, 0.001);

                  for (int ipad = 0; ipad < nPadsR; ++ipad)
                  {
                    const int i = ipad;

                    cER.cd(ipad + 1);
                    gPad->SetLeftMargin(0.12);
                    gPad->SetRightMargin(0.04);
                    gPad->SetBottomMargin(0.12);
                    gPad->SetTopMargin(0.06);

                    if (i < 0 || i >= nPtAll || !perPhoErrRatio[i])
                    {
                      TH1F frame("frame","", 1, 0.0, 2.0);
                      frame.SetMinimum(0.0);
                      frame.SetMaximum(2.0);
                      frame.SetTitle("");
                      frame.GetXaxis()->SetTitle("x_{J}");
                      frame.GetYaxis()->SetTitle("#sigma_{cov} / #sigma_{toy}");
                      frame.Draw("axis");

                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextSize(0.050);
                      tx.DrawLatex(0.16, 0.50, "MISSING");
                      continue;
                    }

                    const PtBin& b = analysisRecoBins[i];

                    TH1* hR = perPhoErrRatio[i];
                    hR->GetXaxis()->SetRangeUser(0.0, 2.0);
                    hR->SetMinimum(0.0);

                    double rMax = 0.0;
                    for (int ib = 1; ib <= hR->GetNbinsX(); ++ib)
                    {
                      const double v = hR->GetBinContent(ib);
                      if (v > rMax) rMax = v;
                    }
                    hR->SetMaximum((rMax > 0.0) ? (1.25 * rMax) : 2.0);
                    hR->SetLineWidth(0);
                    hR->Draw("P E1");

                    // Centered per-pad title
                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(22);
                      tx.SetTextSize(0.042);
                      tx.DrawLatex(0.52, 0.955,
                                   TString::Format("Error ratio, p_{T}^{#gamma} %d-%d GeV, R = %.1f", b.lo, b.hi, R).Data());
                    }

                    // Ensure pad coordinate system is finalized before drawing primitives
                    if (gPad) { gPad->Modified(); gPad->Update(); }

                    // Dashed reference line at y = 1 (persist it in the pad so it appears in saved output)
                    {
                            const double xmin = (gPad ? gPad->GetUxmin() : 0.0);
                            const double xmax = (gPad ? gPad->GetUxmax() : 2.0);
                            TLine* l1 = new TLine(xmin, 1.0, xmax, 1.0);
                            l1->SetLineStyle(2);
                            l1->SetLineWidth(2);
                            l1->SetLineColor(kGray + 2);
                            l1->Draw("same");
                    }

                    // Per-pad cut/trigger label (moved to bottom-right corner)
                    {
                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextAlign(31);
                          tx.SetTextSize(0.038);

                          const double xR = 0.93;
                          tx.DrawLatex(xR, 0.28, TString::Format("Bayes it = %d", kBayesIterXJ).Data());
                          tx.DrawLatex(xR, 0.35, "z_{vtx} < 60 cm");
                          tx.DrawLatex(xR, 0.42, "#Delta #phi > 7#pi/8");
                          tx.DrawLatex(xR, 0.49, "p_{T}^{min, jet} > 5");
                          tx.DrawLatex(xR, 0.56, "Trigger = Photon 4 + MBD NS #geq 1");
                     }
                  }

                  SaveCanvas(cER, JoinPath(toyVsCovOut, "table3x3_ToyUnfoldingVsCovariance_errorRatio.png"));
                }
              }

              bool any = false;
              for (auto* h : perPhoHists) if (h) { any = true; break; }

              if (any)
              {
                  TCanvas c(
                    TString::Format("c_tbl_unf_perPho_%s", rKey.c_str()).Data(),
                    "c_tbl_unf_perPho", 1800, 1300
                  );
                  c.Divide(nPtCols, nPtRows, 0.001, 0.001);

              for (int ipad = 0; ipad < nPtPads; ++ipad)
              {
                const int i = ipad;

                c.cd(ipad + 1);
                gPad->SetLeftMargin(0.12);
                gPad->SetRightMargin(0.04);
                gPad->SetBottomMargin(0.12);
                gPad->SetTopMargin(0.06);

                if (i < 0 || i >= nPtAll)
                {
                  TH1F frame("frame","", 1, 0.0, 2.0);
                  frame.SetMinimum(0.0);
                  frame.SetMaximum(1.0);
                  frame.SetTitle("");
                  frame.GetXaxis()->SetTitle("x_{J}");
                  frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                  frame.Draw("axis");

                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextSize(0.050);
                  tx.DrawLatex(0.16, 0.50, "MISSING");
                  continue;
                }

                const PtBin& b = analysisRecoBins[i];

                if (perPhoHists[i])
                {
                  TH1* h = perPhoHists[i];
                  h->GetXaxis()->SetRangeUser(0.0, 2.0);

                  // Recompute a tight Y-range from actual bin content+error (ignore any previously-set maximum)
                  double maxY = 0.0;
                  const int nxb = h->GetNbinsX();
                  for (int ib = 1; ib <= nxb; ++ib)
                  {
                    const double y  = h->GetBinContent(ib);
                    const double ey = h->GetBinError(ib);
                    const double v  = y + ey;
                    if (v > maxY) maxY = v;
                  }

                  h->SetMinimum(0.0);
                  h->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                  h->Draw("E1");
                }
                else
                {
                  TH1F frame("frame","", 1, 0.0, 2.0);
                  frame.SetMinimum(0.0);
                  frame.SetMaximum(1.0);
                  frame.SetTitle("");
                  frame.GetXaxis()->SetTitle("x_{J}");
                  frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                  frame.Draw("axis");

                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextSize(0.050);
                  tx.DrawLatex(0.16, 0.50, "MISSING");
                }

                // Centered per-pad title
                {
                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextAlign(22);
                  tx.SetTextSize(0.042);

                  tx.DrawLatex(0.52, 0.955,
                               TString::Format("Per-photon particle-level x_{J#gamma}, p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                               b.lo, b.hi, R).Data());
                }

                  // Per-pad cut/trigger label (top-right, 3 lines)
                  {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(31);
                    tx.SetTextSize(0.04);

                      const double xR = 0.93;
                      tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                      tx.DrawLatex(xR, 0.60, TString::Format("Bayes it = %d", kBayesIterXJ).Data());
                      tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                      tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                      tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                  }
              }

                  SaveCanvas(c, JoinPath(finalUnfoldedOut, "table3x3_unfolded_perPhoton_dNdXJ.png"));

                  if (gApplyPurityCorrectionForUnfolding)
                  {
                    bool anyJetEffCorr = false;
                    for (int i = 0; i < nPtAll; ++i)
                    {
                      if (perPhoHists[i] && perPhoHists_jetEffCorr[i])
                      {
                        anyJetEffCorr = true;
                        break;
                      }
                    }

                    if (anyJetEffCorr)
                    {
                      TCanvas cJE(
                        TString::Format("c_tbl_withAndWithoutJetEff_%s", rKey.c_str()).Data(),
                        "c_tbl_withAndWithoutJetEff", 1800, 1300
                      );
                      cJE.Divide(nPtCols, nPtRows, 0.001, 0.001);

                      std::vector<TLegend*> keepLegJE;
                      keepLegJE.reserve((std::size_t)nPtPads);

                      for (int ipad = 0; ipad < nPtPads; ++ipad)
                      {
                        const int i = ipad;

                        cJE.cd(ipad + 1);
                        gPad->SetLeftMargin(0.12);
                        gPad->SetRightMargin(0.04);
                        gPad->SetBottomMargin(0.12);
                        gPad->SetTopMargin(0.06);

                        if (i < 0 || i >= nPtAll)
                        {
                          TH1F frame("frame","", 1, 0.0, 2.0);
                          frame.SetMinimum(0.0);
                          frame.SetMaximum(1.0);
                          frame.SetTitle("");
                          frame.GetXaxis()->SetTitle("x_{J}");
                          frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                          frame.Draw("axis");

                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextSize(0.050);
                          tx.DrawLatex(0.16, 0.50, "MISSING");
                          continue;
                        }

                        const PtBin& b = analysisRecoBins[i];

                        TH1* hNo  = perPhoHists[i];
                        TH1* hYes = perPhoHists_jetEffCorr[i];

                        if (hNo && hYes)
                        {
                          hNo->GetXaxis()->SetRangeUser(0.0, 2.0);
                          hYes->GetXaxis()->SetRangeUser(0.0, 2.0);

                          hNo->SetMarkerStyle(20);
                          hNo->SetMarkerColor(kBlack);
                          hNo->SetLineColor(kBlack);
                          hNo->SetLineWidth(2);

                          hYes->SetMarkerStyle(24);
                          hYes->SetMarkerColor(kRed + 1);
                          hYes->SetLineColor(kRed + 1);
                          hYes->SetLineWidth(2);

                          double maxY = 0.0;
                          const int nxb = hNo->GetNbinsX();
                          for (int ib = 1; ib <= nxb; ++ib)
                          {
                            const double v1 = hNo->GetBinContent(ib)  + hNo->GetBinError(ib);
                            const double v2 = hYes->GetBinContent(ib) + hYes->GetBinError(ib);
                            if (v1 > maxY) maxY = v1;
                            if (v2 > maxY) maxY = v2;
                          }

                          hNo->SetMinimum(0.0);
                          hNo->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                          hNo->Draw("E1");
                          hYes->Draw("E1 same");

                          TLegend* leg = new TLegend(0.54, 0.33, 0.90, 0.55);
                          leg->SetBorderSize(0);
                          leg->SetFillStyle(0);
                          leg->SetTextFont(42);
                          leg->SetTextSize(0.038);
                          leg->AddEntry(hNo,  "unfolded data (no jet eff. corr.)", "pe");
                          leg->AddEntry(hYes, "unfolded data (+ jet eff. corr.)",  "pe");
                          leg->Draw();
                          keepLegJE.push_back(leg);
                        }
                        else
                        {
                          TH1F frame("frame","", 1, 0.0, 2.0);
                          frame.SetMinimum(0.0);
                          frame.SetMaximum(1.0);
                          frame.SetTitle("");
                          frame.GetXaxis()->SetTitle("x_{J}");
                          frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                          frame.Draw("axis");

                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextSize(0.050);
                          tx.DrawLatex(0.16, 0.50, "MISSING");
                        }

                        {
                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextAlign(22);
                          tx.SetTextSize(0.040);

                          tx.DrawLatex(0.52, 0.955,
                                       TString::Format("Without vs with jet eff. corr., p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                                       b.lo, b.hi, R).Data());
                        }

                        {
                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextAlign(31);
                          tx.SetTextSize(0.04);

                          const double xR = 0.93;
                          tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                          tx.DrawLatex(xR, 0.60, TString::Format("Bayes it = %d", kBayesIterXJ).Data());
                          tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                          tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                          tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                        }
                      }

                      SaveCanvas(cJE, JoinPath(withAndWithoutJetEffOut, "table3x3_withAndWithoutJetEffCorr.png"));

                      for (auto* p : keepLegJE) { delete p; }

                      for (int i = 0; i < nPtAll; ++i)
                      {
                        const PtBin& b = analysisRecoBins[i];

                        TH1* hNo  = perPhoHists[i];
                        TH1* hYes = perPhoHists_jetEffCorr[i];
                        if (!hNo || !hYes) continue;

                        TCanvas cSingle(
                          TString::Format("c_withAndWithoutJetEff_%s_%d", rKey.c_str(), i + 1).Data(),
                          "c_withAndWithoutJetEff", 900, 700
                        );
                        ApplyCanvasMargins1D(cSingle);

                        hNo->GetXaxis()->SetRangeUser(0.0, 2.0);
                        hYes->GetXaxis()->SetRangeUser(0.0, 2.0);

                        hNo->SetMarkerStyle(20);
                        hNo->SetMarkerColor(kBlack);
                        hNo->SetLineColor(kBlack);
                        hNo->SetLineWidth(2);

                        hYes->SetMarkerStyle(24);
                        hYes->SetMarkerColor(kRed + 1);
                        hYes->SetLineColor(kRed + 1);
                        hYes->SetLineWidth(2);

                        double maxY = 0.0;
                        const int nxb = hNo->GetNbinsX();
                        for (int ib = 1; ib <= nxb; ++ib)
                        {
                          const double v1 = hNo->GetBinContent(ib)  + hNo->GetBinError(ib);
                          const double v2 = hYes->GetBinContent(ib) + hYes->GetBinError(ib);
                          if (v1 > maxY) maxY = v1;
                          if (v2 > maxY) maxY = v2;
                        }

                        hNo->SetMinimum(0.0);
                        hNo->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                        hNo->Draw("E1");
                        hYes->Draw("E1 same");

                        TLegend leg(0.58, 0.35, 0.90, 0.53);
                        leg.SetBorderSize(0);
                        leg.SetFillStyle(0);
                        leg.SetTextFont(42);
                        leg.SetTextSize(0.03);
                        leg.AddEntry(hNo,  "unfolded data (no jet eff. corr.)", "pe");
                        leg.AddEntry(hYes, "unfolded data (+ jet eff. corr.)",  "pe");
                        leg.Draw();

                        {
                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextAlign(22);
                          tx.SetTextSize(0.040);

                          tx.DrawLatex(0.50, 0.955,
                                       TString::Format("Without vs with jet eff. corr., p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                                       b.lo, b.hi, R).Data());
                        }

                        {
                          TLatex tx;
                          tx.SetNDC();
                          tx.SetTextFont(42);
                          tx.SetTextAlign(31);
                          tx.SetTextSize(0.04);

                          const double xR = 0.93;
                          tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                          tx.DrawLatex(xR, 0.60, TString::Format("Bayes it = %d", kBayesIterXJ).Data());
                          tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                          tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                          tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                        }

                        SaveCanvas(
                          cSingle,
                          JoinPath(
                            withAndWithoutJetEffOut,
                            TString::Format("xJ_withAndWithoutJetEffCorr_pTbin%d.png", i + 1).Data()
                          )
                        );
                      }
                    }
                  }

                  // -------------------------------------------------------------------
                  // 2x4 summary table: before vs after unfolding (DATA)
                  //   output: <rOut>/before_after_unfoldingOverlay_data/table2x4_before_after_unfoldingOverlay_data.png
                  // -------------------------------------------------------------------
              {
                  TCanvas cBA(
                    TString::Format("c_tbl_beforeAfter_data_%s", rKey.c_str()).Data(),
                    "c_tbl_beforeAfter_data", 1800, 1300
                  );
                  cBA.Divide(nPtCols, nPtRows, 0.001, 0.001);

                  std::vector<TLegend*> keepLegBA;
                  keepLegBA.reserve((std::size_t)nPtPads);

                  for (int ipad = 0; ipad < nPtPads; ++ipad)
                  {
                    const int i = ipad;

                    cBA.cd(ipad + 1);
                    gPad->SetLeftMargin(0.12);
                    gPad->SetRightMargin(0.04);
                    gPad->SetBottomMargin(0.12);
                    gPad->SetTopMargin(0.06);

                    if (i < 0 || i >= nPtAll)
                    {
                      TH1F frame("frame","", 1, 0.0, 2.0);
                      frame.SetMinimum(0.0);
                      frame.SetMaximum(1.0);
                      frame.SetTitle("");
                      frame.GetXaxis()->SetTitle("x_{J}");
                      frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                      frame.Draw("axis");

                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextSize(0.050);
                      tx.DrawLatex(0.16, 0.50, "MISSING");
                      continue;
                    }

                    const PtBin& b = analysisRecoBins[i];

                    TH1* hA = perPhoBeforeDataHists[i];
                    TH1* hB = perPhoHists[i];

                    if (hA && hB)
                    {
                        hA->GetXaxis()->SetRangeUser(0.0, 2.0);
                        hB->GetXaxis()->SetRangeUser(0.0, 2.0);

                        hA->SetMarkerStyle(24);
                        hA->SetMarkerColor(kRed + 1);
                        hA->SetLineColor(kRed + 1);
                        hA->SetLineWidth(2);

                        hB->SetMarkerStyle(20);
                        hB->SetMarkerColor(kBlue + 1);
                        hB->SetLineColor(kBlue + 1);
                        hB->SetLineWidth(2);

                        double maxY = 0.0;
                        const int nxb = hA->GetNbinsX();
                        for (int ib = 1; ib <= nxb; ++ib)
                        {
                          const double v1 = hA->GetBinContent(ib) + hA->GetBinError(ib);
                          const double v2 = hB->GetBinContent(ib) + hB->GetBinError(ib);
                          if (v1 > maxY) maxY = v1;
                          if (v2 > maxY) maxY = v2;
                        }

                        hA->SetMinimum(0.0);
                        hA->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                        hA->Draw("E1");
                        hB->Draw("E1 same");

                        TLegend* leg = new TLegend(0.65, 0.4, 0.88, 0.55);
                        leg->SetBorderSize(0);
                        leg->SetFillStyle(0);
                        leg->SetTextFont(42);
                        leg->SetTextSize(0.040);
                        leg->AddEntry(hA, "before unfolding data", "pe");
                        leg->AddEntry(hB, "unfolded data",        "pe");
                        leg->Draw();
                        keepLegBA.push_back(leg);
                    }
                    else
                    {
                      TH1F frame("frame","", 1, 0.0, 2.0);
                      frame.SetMinimum(0.0);
                      frame.SetMaximum(1.0);
                      frame.SetTitle("");
                      frame.GetXaxis()->SetTitle("x_{J}");
                      frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                      frame.Draw("axis");

                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextSize(0.050);
                      tx.DrawLatex(0.16, 0.50, "MISSING");
                    }

                    // Centered per-pad title
                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(22);
                      tx.SetTextSize(0.042);

                      tx.DrawLatex(0.52, 0.955,
                                   TString::Format("Before vs after unfolding (DATA), p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                                   b.lo, b.hi, R).Data());
                    }

                    // Per-pad cut/trigger label (top-right, 3 lines)
                    {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(31);
                      tx.SetTextSize(0.04);

                        const double xR = 0.93;
                        tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                        tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                        tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                        tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                    }
                  }

                  SaveCanvas(cBA, JoinPath(beforeAfterDataOut, "table3x3_before_after_unfoldingOverlay_data.png"));

                  for (auto* p : keepLegBA) { delete p; }
                }

                // -------------------------------------------------------------------
                // NEW: 2x4 summary table: truth MC vs unfolded data
                //   output: <rOut>/before_after_unfoldingOverlay_truth/table2x4_before_after_unfoldingOverlay_truth.png
                // -------------------------------------------------------------------
                {
                  TCanvas cTU(
                      TString::Format("c_tbl_truthVsUnf_%s", rKey.c_str()).Data(),
                      "c_tbl_truthVsUnf", 1800, 1300
                  );
                  cTU.Divide(nPtCols, nPtRows, 0.001, 0.001);

                  std::vector<TLegend*> keepLegTU;
                  keepLegTU.reserve((std::size_t)nPtPads);

                  for (int ipad = 0; ipad < nPtPads; ++ipad)
                  {
                      const int i = ipad;

                      cTU.cd(ipad + 1);
                      gPad->SetLeftMargin(0.12);
                      gPad->SetRightMargin(0.04);
                      gPad->SetBottomMargin(0.12);
                      gPad->SetTopMargin(0.06);

                      if (i < 0 || i >= nPtAll)
                      {
                        TH1F frame("frame","", 1, 0.0, 2.0);
                        frame.SetMinimum(0.0);
                        frame.SetMaximum(1.0);
                        frame.SetTitle("");
                        frame.GetXaxis()->SetTitle("x_{J}");
                        frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                        frame.Draw("axis");

                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextSize(0.050);
                        tx.DrawLatex(0.16, 0.50, "MISSING");
                        continue;
                      }

                      const PtBin& b = analysisRecoBins[i];

                      TH1* hT = perPhoTruthHists[i];
                      TH1* hU = perPhoHists[i];

                      if (hT && hU)
                      {
                        hT->GetXaxis()->SetRangeUser(0.0, 2.0);
                        hU->GetXaxis()->SetRangeUser(0.0, 2.0);

                        double maxY = 0.0;
                        const int nxb = hT->GetNbinsX();
                        for (int ib = 1; ib <= nxb; ++ib)
                        {
                          const double v1 = hT->GetBinContent(ib) + hT->GetBinError(ib);
                          const double v2 = hU->GetBinContent(ib) + hU->GetBinError(ib);
                          if (v1 > maxY) maxY = v1;
                          if (v2 > maxY) maxY = v2;
                        }

                        hT->SetMinimum(0.0);
                        hT->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                        hT->Draw("E1");
                        hU->SetMarkerStyle(20);
                        hU->SetMarkerColor(kBlue);
                        hU->SetLineColor(kBlue);
                        hU->Draw("E1 same");

                        TLegend* leg = new TLegend(0.69, 0.33, 0.91, 0.58);
                        leg->SetBorderSize(0);
                        leg->SetFillStyle(0);
                        leg->SetTextFont(42);
                        leg->SetTextSize(0.038);
                        leg->AddEntry(hT, "truth MC",      "pe");
                        leg->AddEntry(hU, "unfolded data", "pe");
                        leg->Draw();
                        keepLegTU.push_back(leg);
                      }
                      else
                      {
                        TH1F frame("frame","", 1, 0.0, 2.0);
                        frame.SetMinimum(0.0);
                        frame.SetMaximum(1.0);
                        frame.SetTitle("");
                        frame.GetXaxis()->SetTitle("x_{J}");
                        frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                        frame.Draw("axis");

                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextSize(0.050);
                        tx.DrawLatex(0.16, 0.50, "MISSING");
                      }

                      // Centered per-pad title
                      {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(22);
                        tx.SetTextSize(0.042);

                        tx.DrawLatex(0.52, 0.955,
                                     TString::Format("Truth MC vs unfolded data, p_{T}^{#gamma} %d-%d GeV, R = %.1f",
                                                     b.lo, b.hi, R).Data());
                      }

                      // Per-pad cut/trigger label (top-right, 3 lines)
                      {
                        TLatex tx;
                        tx.SetNDC();
                        tx.SetTextFont(42);
                        tx.SetTextAlign(31);
                        tx.SetTextSize(0.04);

                          const double xR = 0.93;
                          tx.DrawLatex(xR, 0.67, "z_{vtx} < 60 cm");
                          tx.DrawLatex(xR, 0.74, "#Delta #phi > 7#pi/8");
                          tx.DrawLatex(xR, 0.81, "p_{T}^{min, jet} > 5");
                          tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                      }
                  }

                  SaveCanvas(cTU, JoinPath(beforeAfterTruthOut, "table3x3_before_after_unfoldingOverlay_truth.png"));

                  for (auto* p : keepLegTU) { delete p; }
               }
            }
          }

          // ----------------------------------------------------------------------
          // Iteration stability / regularization diagnostic (DATA):
          //   • total relative deviation between successive iterations
          //   • total relative statistical uncertainty of the unfolded spectrum
          //
          // Output:
          //   <rOut>/unfold_iterStability_relChange_relStat.png
          //   <rOut>/unfold_iterStability_quadratureSum.png
          //   <rOut>/unfold_iterStability_bestIteration.txt
          // ----------------------------------------------------------------------
          {
            const int kMaxIterPlot = kMaxBayesIterXJScan;

            cout << ANSI_BOLD_CYN
                 << "[UNF ITER QA] Building iteration-stability plot from cached scan\n"
                 << "  rKey=" << rKey << "  R=" << std::fixed << std::setprecision(1) << R << "\n"
                 << "  outdir=" << rOut << "\n"
                 << "  kMaxIterPlot=" << kMaxIterPlot << "\n"
                 << "  selectedBayesIter=" << kBayesIterXJ << "\n"
                 << ANSI_RESET;

            const std::vector<double>& xIt         = xjIterScan.xIt;
            const std::vector<double>& exIt        = xjIterScan.exIt;
            const std::vector<double>& yRelStat    = xjIterScan.yRelStat;
            const std::vector<double>& eyRelStat   = xjIterScan.eyRelStat;
            const std::vector<double>& yRelChange  = xjIterScan.yRelChange;
            const std::vector<double>& eyRelChange = xjIterScan.eyRelChange;
            const std::vector<double>& yQuad       = xjIterScan.yQuad;
            const std::vector<double>& eyQuad      = xjIterScan.eyQuad;
            const int bestIt                       = xjIterScan.bestIt;
            const double bestQuad                  = xjIterScan.bestQuad;

            cout << ANSI_BOLD_CYN
                 << "[UNF ITER QA] Points prepared: n=" << xIt.size()
                 << " (expected " << kMaxIterPlot << " for it=1..kMax)\n"
                 << ANSI_RESET;

              if (!xIt.empty())
              {
                vector<double> xItNo1;
                vector<double> exItNo1;
                vector<double> yRelStatNo1;
                vector<double> eyRelStatNo1;
                vector<double> yRelChangeNo1;
                vector<double> eyRelChangeNo1;
                vector<double> yQuadNo1;
                vector<double> eyQuadNo1;

                for (size_t i = 0; i < xIt.size(); ++i)
                {
                  if (xIt[i] <= 1.0) continue;
                  xItNo1.push_back(xIt[i]);
                  exItNo1.push_back(exIt[i]);
                  yRelStatNo1.push_back(yRelStat[i]);
                  eyRelStatNo1.push_back(eyRelStat[i]);
                  yRelChangeNo1.push_back(yRelChange[i]);
                  eyRelChangeNo1.push_back(eyRelChange[i]);
                  yQuadNo1.push_back(yQuad[i]);
                  eyQuadNo1.push_back(eyQuad[i]);
                }

                double yMax = 0.0;
                for (size_t i = 0; i < yRelStat.size();   ++i) yMax = std::max(yMax, yRelStat[i]);
                for (size_t i = 0; i < yRelChange.size(); ++i) yMax = std::max(yMax, yRelChange[i]);
                if (yMax <= 0.0) yMax = 1.0;

                cout << ANSI_BOLD_CYN << "[UNF ITER QA] yMax=" << std::setprecision(6) << yMax << ANSI_RESET << "\n";

                TCanvas cIt(TString::Format("c_iterStability_%s", rKey.c_str()).Data(), "c_iterStability", 900, 700);
                ApplyCanvasMargins1D(cIt);

                TH1F frame("frame","", 1, 1.0, (double)kMaxIterPlot + 0.5);
                frame.SetMinimum(0.0);
                frame.SetMaximum(1.20 * yMax);
                frame.SetTitle("");
                frame.GetXaxis()->SetTitle("Iteration");
                frame.GetYaxis()->SetTitle("Relative quantity");
                frame.Draw("axis");

                TGraphErrors gStat((int)xIt.size(), &xIt[0], &yRelStat[0], &exIt[0], &eyRelStat[0]);
                gStat.SetMarkerStyle(20);
                gStat.SetMarkerSize(1.1);
                gStat.SetMarkerColor(kRed + 1);
                gStat.SetLineColor(kRed + 1);
                gStat.SetLineWidth(2);
                gStat.Draw("P same");

                TGraphErrors gChg((int)xIt.size(), &xIt[0], &yRelChange[0], &exIt[0], &eyRelChange[0]);
                gChg.SetMarkerStyle(20);
                gChg.SetMarkerSize(1.1);
                gChg.SetMarkerColor(kBlue);
                gChg.SetLineColor(kBlue);
                gChg.SetLineWidth(2);
                gChg.Draw("P same");

                TLegend leg(0.14, 0.78, 0.54, 0.92);
                leg.SetBorderSize(0);
                leg.SetFillStyle(0);
                leg.SetTextFont(42);
                leg.SetTextSize(0.032);
                leg.AddEntry(&gStat, "total relative stat. uncertainty", "p");
                leg.AddEntry(&gChg,  "total relative deviation (it vs it-1)", "p");
                leg.Draw();

                {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(13);
                    tx.SetTextSize(0.032);
                    tx.DrawLatex(0.14, 0.74, "2D (p_{T}^{#gamma}, x_{J}) unfolding");
                }

                {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(22);
                    tx.SetTextSize(0.040);
                    tx.DrawLatex(0.50, 0.965, TString::Format("Iteration Stability, R = %.1f, Photon 4 + MBD NS #geq 1, Run24pp", R).Data());
                }

                cIt.Modified();
                cIt.Update();

                const std::string outPng = JoinPath(iterationStabilityOut, "unfold_iterStability_relChange_relStat.png");
                cout << ANSI_BOLD_CYN << "[UNF ITER QA] Saving: " << outPng << ANSI_RESET << "\n";
                SaveCanvas(cIt, outPng);

                if (!xItNo1.empty())
                {
                  double yMaxNo1 = 0.0;
                  for (size_t i = 0; i < yRelStatNo1.size();   ++i) yMaxNo1 = std::max(yMaxNo1, yRelStatNo1[i]);
                  for (size_t i = 0; i < yRelChangeNo1.size(); ++i) yMaxNo1 = std::max(yMaxNo1, yRelChangeNo1[i]);
                  if (yMaxNo1 <= 0.0) yMaxNo1 = 1.0;

                  TCanvas cItNo1(TString::Format("c_iterStability_noIter1_%s", rKey.c_str()).Data(), "c_iterStability_noIter1", 900, 700);
                  ApplyCanvasMargins1D(cItNo1);

                  TH1F frameNo1("frameNo1","", 1, 2.0, (double)kMaxIterPlot + 0.5);
                  frameNo1.SetMinimum(0.0);
                  frameNo1.SetMaximum(1.20 * yMaxNo1);
                  frameNo1.SetTitle("");
                  frameNo1.GetXaxis()->SetTitle("Iteration");
                  frameNo1.GetYaxis()->SetTitle("Relative quantity");
                  frameNo1.Draw("axis");

                  TGraphErrors gStatNo1((int)xItNo1.size(), &xItNo1[0], &yRelStatNo1[0], &exItNo1[0], &eyRelStatNo1[0]);
                  gStatNo1.SetMarkerStyle(20);
                  gStatNo1.SetMarkerSize(1.1);
                  gStatNo1.SetMarkerColor(kRed + 1);
                  gStatNo1.SetLineColor(kRed + 1);
                  gStatNo1.SetLineWidth(2);
                  gStatNo1.Draw("P same");

                  TGraphErrors gChgNo1((int)xItNo1.size(), &xItNo1[0], &yRelChangeNo1[0], &exItNo1[0], &eyRelChangeNo1[0]);
                  gChgNo1.SetMarkerStyle(20);
                  gChgNo1.SetMarkerSize(1.1);
                  gChgNo1.SetMarkerColor(kBlue);
                  gChgNo1.SetLineColor(kBlue);
                  gChgNo1.SetLineWidth(2);
                  gChgNo1.Draw("P same");

                  TLegend legNo1(0.14, 0.78, 0.54, 0.92);
                  legNo1.SetBorderSize(0);
                  legNo1.SetFillStyle(0);
                  legNo1.SetTextFont(42);
                  legNo1.SetTextSize(0.032);
                  legNo1.AddEntry(&gStatNo1, "total relative stat. uncertainty", "p");
                  legNo1.AddEntry(&gChgNo1,  "total relative deviation (it vs it-1)", "p");
                  legNo1.Draw();

                  {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(13);
                      tx.SetTextSize(0.032);
                      tx.DrawLatex(0.14, 0.74, "2D (p_{T}^{#gamma}, x_{J}) unfolding (it #geq 2 only)");
                  }

                  {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(22);
                      tx.SetTextSize(0.040);
                      tx.DrawLatex(0.50, 0.965, TString::Format("Iteration Stability, R = %.1f, Photon 4 + MBD NS #geq 1, Run24pp", R).Data());
                  }

                  cItNo1.Modified();
                  cItNo1.Update();

                  const std::string outPngNo1 = JoinPath(iterationStabilityOut, "unfold_iterStability_relChange_relStat_noIter1.png");
                  cout << ANSI_BOLD_CYN << "[UNF ITER QA] Saving: " << outPngNo1 << ANSI_RESET << "\n";
                  SaveCanvas(cItNo1, outPngNo1);
                }

                double yMaxQ = 0.0;
                for (size_t i = 0; i < yQuad.size(); ++i) yMaxQ = std::max(yMaxQ, yQuad[i]);
                if (yMaxQ <= 0.0) yMaxQ = 0.1;

                TCanvas cQuad(TString::Format("c_iterQuadrature_%s", rKey.c_str()).Data(), "c_iterQuadrature", 900, 700);
                ApplyCanvasMargins1D(cQuad);

                TH1F frameQ("frameQ","", 1, 1.0, (double)kMaxIterPlot + 0.5);
                frameQ.SetMinimum(0.0);
                frameQ.SetMaximum(1.20 * yMaxQ);
                frameQ.SetTitle("");
                frameQ.GetXaxis()->SetTitle("Iteration");
                frameQ.GetYaxis()->SetTitle("Quadrature sum");
                frameQ.Draw("axis");

                TGraphErrors gQuad((int)xIt.size(), &xIt[0], &yQuad[0], &exIt[0], &eyQuad[0]);
                gQuad.SetMarkerStyle(20);
                gQuad.SetMarkerSize(1.1);
                gQuad.SetMarkerColor(kBlack);
                gQuad.SetLineColor(kBlack);
                gQuad.SetLineWidth(2);
                gQuad.Draw("LP same");

                TGraphErrors gBest;
                if (bestIt > 0 && std::isfinite(bestQuad))
                {
                  const double bestX[1]  = { (double)bestIt };
                  const double bestY[1]  = { bestQuad };
                  const double bestEX[1] = { 0.0 };
                  const double bestEY[1] = { 0.0 };

                  gBest = TGraphErrors(1, bestX, bestY, bestEX, bestEY);
                  gBest.SetMarkerStyle(29);
                  gBest.SetMarkerSize(1.6);
                  gBest.SetMarkerColor(kRed + 1);
                  gBest.SetLineColor(kRed + 1);
                  gBest.Draw("P same");
                }

                TLegend legQ(0.14, 0.78, 0.54, 0.92);
                legQ.SetBorderSize(0);
                legQ.SetFillStyle(0);
                legQ.SetTextFont(42);
                legQ.SetTextSize(0.032);
                legQ.AddEntry(&gQuad, "quadrature sum", "lp");
                if (bestIt > 0 && std::isfinite(bestQuad))
                {
                  legQ.AddEntry(&gBest, TString::Format("minimum: it=%d", bestIt).Data(), "p");
                }
                legQ.Draw();

                {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(13);
                    tx.SetTextSize(0.032);
                    tx.DrawLatex(0.14, 0.74, "2D (p_{T}^{#gamma}, x_{J}) unfolding");
                }

                {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(22);
                    tx.SetTextSize(0.040);
                    tx.DrawLatex(0.50, 0.965, TString::Format("Quadrature-sum optimization, R = %.1f, Photon 4 + MBD NS #geq 1, Run24pp", R).Data());
                }

                const std::string outQuadPng = JoinPath(iterationStabilityOut, "unfold_iterStability_quadratureSum.png");
                cout << ANSI_BOLD_CYN << "[UNF ITER QA] Saving: " << outQuadPng << ANSI_RESET << "\n";
                SaveCanvas(cQuad, outQuadPng);

                if (!xItNo1.empty())
                {
                  double yMaxQNo1 = 0.0;
                  for (size_t i = 0; i < yQuadNo1.size(); ++i) yMaxQNo1 = std::max(yMaxQNo1, yQuadNo1[i]);
                  if (yMaxQNo1 <= 0.0) yMaxQNo1 = 0.1;

                  TCanvas cQuadNo1(TString::Format("c_iterQuadrature_noIter1_%s", rKey.c_str()).Data(), "c_iterQuadrature_noIter1", 900, 700);
                  ApplyCanvasMargins1D(cQuadNo1);

                  TH1F frameQNo1("frameQNo1","", 1, 2.0, (double)kMaxIterPlot + 0.5);
                  frameQNo1.SetMinimum(0.0);
                  frameQNo1.SetMaximum(1.20 * yMaxQNo1);
                  frameQNo1.SetTitle("");
                  frameQNo1.GetXaxis()->SetTitle("Iteration");
                  frameQNo1.GetYaxis()->SetTitle("Quadrature sum");
                  frameQNo1.Draw("axis");

                  TGraphErrors gQuadNo1((int)xItNo1.size(), &xItNo1[0], &yQuadNo1[0], &exItNo1[0], &eyQuadNo1[0]);
                  gQuadNo1.SetMarkerStyle(20);
                  gQuadNo1.SetMarkerSize(1.1);
                  gQuadNo1.SetMarkerColor(kBlack);
                  gQuadNo1.SetLineColor(kBlack);
                  gQuadNo1.SetLineWidth(2);
                  gQuadNo1.Draw("LP same");

                  TLegend legQNo1(0.14, 0.78, 0.54, 0.92);
                  legQNo1.SetBorderSize(0);
                  legQNo1.SetFillStyle(0);
                  legQNo1.SetTextFont(42);
                  legQNo1.SetTextSize(0.032);
                  legQNo1.AddEntry(&gQuadNo1, "quadrature sum", "lp");
                  legQNo1.Draw();

                  {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(13);
                      tx.SetTextSize(0.032);
                      tx.DrawLatex(0.14, 0.74, "2D (p_{T}^{#gamma}, x_{J}) unfolding (it #geq 2 only)");
                  }

                  {
                      TLatex tx;
                      tx.SetNDC();
                      tx.SetTextFont(42);
                      tx.SetTextAlign(22);
                      tx.SetTextSize(0.040);
                      tx.DrawLatex(0.50, 0.965, TString::Format("Quadrature-sum optimization, R = %.1f, Photon 4 + MBD NS #geq 1, Run24pp", R).Data());
                  }

                  const std::string outQuadPngNo1 = JoinPath(iterationStabilityOut, "unfold_iterStability_quadratureSum_noIter1.png");
                  cout << ANSI_BOLD_CYN << "[UNF ITER QA] Saving: " << outQuadPngNo1 << ANSI_RESET << "\n";
                  SaveCanvas(cQuadNo1, outQuadPngNo1);
                }

                vector<string> iterSummary;
                iterSummary.push_back("xJ iteration-stability summary");
                iterSummary.push_back(TString::Format("radius = %s (R=%.1f)", rKey.c_str(), R).Data());
                iterSummary.push_back(TString::Format("best iteration from quadrature sum = %d", bestIt).Data());
                iterSummary.push_back(TString::Format("minimum quadrature sum = %.10g", bestQuad).Data());
                iterSummary.push_back("relChange is defined as the successive-iterate difference it vs (it-1), with iteration 1 using the explicit 0->1 baseline.");
                iterSummary.push_back("iteration 0 baseline is the measured reco global-bin histogram mapped onto the truth comparison axis, so iteration 1 is included in the best-iteration choice.");
                iterSummary.push_back("additional plots with suffix _noIter1 exclude iteration 1 from the display only; they do not change the scan or best-iteration choice.");
                iterSummary.push_back("The stability scan is evaluated only in the 9 analysis p_{T}^{#gamma} bins and in-range x_{J} bins; support / underflow / overflow bins are excluded.");
                iterSummary.push_back("relStat is built from the full RooUnfold::kCovToy covariance restricted to those physics bins.");
                iterSummary.push_back("quadrature sum definition: sqrt(relStat^2 + relChange^2)");
                iterSummary.push_back("");

                for (size_t i = 0; i < xIt.size(); ++i)
                {
                  iterSummary.push_back(
                    TString::Format("it=%d  relStat=%.10g  relChange=%.10g  quadratureSum=%.10g",
                                    (int)xIt[i], yRelStat[i], yRelChange[i], yQuad[i]).Data()
                  );
                }

                WriteTextFile(JoinPath(rOut, "unfold_iterStability_bestIteration.txt"), iterSummary);

                cout << ANSI_BOLD_CYN
                     << "[UNF ITER QA] best iteration from quadrature sum = " << bestIt
                     << "  minimum = " << bestQuad
                     << ANSI_RESET << "\n";
              }
              else
              {
                cout << ANSI_BOLD_RED << "[UNF ITER QA][WARN] No points to plot; not saving PNG." << ANSI_RESET << "\n";
              }
          }

          WriteTextFile(JoinPath(rOut, "summary_rooUnfold_pipeline.txt"), lines);

          {
            const string outRoot = JoinPath(rOut, "rooUnfold_outputs.root");
            TFile f(outRoot.c_str(), "RECREATE");
            if (f.IsOpen())
            {
                if (hPhoUnfoldTruth)      hPhoUnfoldTruth->Write("h_phoTruth_unfolded_data");
                if (hPhoUnfoldTruth_cov)  hPhoUnfoldTruth_cov->Write("h_phoTruth_unfolded_data_covariance");
                if (h2RspSim)             h2RspSim->Write("h2_rsp_truthVsReco_global");
                if (hRsp_measXtruth)      hRsp_measXtruth->Write("h2_rsp_recoVsTruth_global");
                if (hMeasDataGlob)        hMeasDataGlob->Write("h_measData_global");
                if (hMeasSimGlob)         hMeasSimGlob->Write("h_measSim_global");
                if (hTruthSimGlob)        hTruthSimGlob->Write("h_truthSim_global");
                if (hUnfoldTruthGlob)     hUnfoldTruthGlob->Write("h_truthUnfold_global");
                if (hUnfoldTruthGlob_cov) hUnfoldTruthGlob_cov->Write("h_truthUnfold_global_covariance");
                if (h2UnfoldTruth)             h2UnfoldTruth->Write("h2_truthUnfold_pTgamma_xJ");
                if (h2UnfoldTruth_cov)         h2UnfoldTruth_cov->Write("h2_truthUnfold_pTgamma_xJ_covariance");
                if (h2UnfoldTruth_jetEffCorr)  h2UnfoldTruth_jetEffCorr->Write("h2_truthUnfold_pTgamma_xJ_jetEffCorr");
                if (h2JetEff)                  h2JetEff->Write("h2_jetEff_pTgamma_xJ");

                for (int i = 0; i < nPtAll; ++i)
                {
                    if (perPhoHists[i])             perPhoHists[i]->Write(TString::Format("h_xJ_unf_perPho_pTbin%d", i + 1).Data());
                    if (perPhoHists_cov[i])         perPhoHists_cov[i]->Write(TString::Format("h_xJ_unf_perPho_cov_pTbin%d", i + 1).Data());
                    if (perPhoHists_jetEffCorr[i])  perPhoHists_jetEffCorr[i]->Write(TString::Format("h_xJ_unf_perPho_jetEffCorr_pTbin%d", i + 1).Data());
                }

              f.Close();
            }
          }

        // Keep clones for an all-radii overlay table produced after the per-radius loop
        {
            auto& vv = perPhoHistsByRKey[rKey];
            vv.assign(nPtAll, nullptr);

            for (int i = 0; i < nPtAll; ++i)
            {
              if (!perPhoHists[i]) { vv[i] = nullptr; continue; }

              vv[i] = (TH1*)perPhoHists[i]->Clone(
                TString::Format("%s_clone_radiiOverlay_pTbin%d", perPhoHists[i]->GetName(), i + 1).Data()
              );
              if (vv[i])
              {
                vv[i]->SetDirectory(nullptr);
                EnsureSumw2(vv[i]);
              }
            }
          }

          auto DrawRatioOverlaySummary =
            [&](const vector<TH1*>& hs,
                const string& outDir,
                const string& outName,
                const string& yTitle,
                const string& titleText)->void
          {
            bool anyRatio = false;
            for (auto* h : hs) if (h) { anyRatio = true; break; }
            if (!anyRatio) return;

            vector<int> colors =
            {
              kBlack,
              kRed + 1,
              kBlue + 1,
              kGreen + 2,
              kMagenta + 1,
              kOrange + 7,
              kCyan + 2,
              kViolet + 1,
              kPink + 7
            };

            double maxY = 0.0;
            for (int i = 0; i < nPtAll; ++i)
            {
              if (!hs[i]) continue;
              hs[i]->GetXaxis()->SetRangeUser(0.0, 2.0);
              for (int ib = 1; ib <= hs[i]->GetNbinsX(); ++ib)
              {
                const double v = hs[i]->GetBinContent(ib) + hs[i]->GetBinError(ib);
                if (std::isfinite(v) && v > maxY) maxY = v;
              }
            }

            TCanvas c(TString::Format("c_ratioOverlay_%s_%s", rKey.c_str(), outName.c_str()).Data(),
                      "c_ratioOverlay", 950, 750);
            ApplyCanvasMargins1D(c);

            TH1F frame("frame_ratio_overlay","", 1, 0.0, 2.0);
            frame.SetMinimum(0.0);
            frame.SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 2.0);
            frame.SetTitle("");
            frame.GetXaxis()->SetTitle("x_{J}");
            frame.GetYaxis()->SetTitle(yTitle.c_str());
            frame.Draw("axis");

            TLine l1(0.0, 1.0, 2.0, 1.0);
            l1.SetLineStyle(2);
            l1.SetLineWidth(2);
            l1.Draw("same");

            TLegend leg(0.58, 0.52, 0.90, 0.88);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.SetTextFont(42);
            leg.SetTextSize(0.032);

            for (int i = 0; i < nPtAll; ++i)
            {
              if (!hs[i]) continue;

              TH1* h = hs[i];
              h->SetLineColor(colors[(std::size_t)i % colors.size()]);
              h->SetMarkerColor(colors[(std::size_t)i % colors.size()]);
              h->SetMarkerStyle(20);
              h->SetMarkerSize(0.95);
              h->SetLineWidth(2);
              h->Draw("E1 X0 same");

              const PtBin& b = analysisRecoBins[i];
              leg.AddEntry(h, TString::Format("%d-%d GeV", b.lo, b.hi).Data(), "pe");
            }

            leg.Draw();

            TLatex tx;
            tx.SetNDC();
            tx.SetTextFont(42);
            tx.SetTextAlign(22);
            tx.SetTextSize(0.040);
            tx.DrawLatex(0.50, 0.965, titleText.c_str());

            SaveCanvas(c, JoinPath(outDir, outName));
          };

          DrawRatioOverlaySummary(
            ratioBeforeVsAfterHists,
            beforeAfterDataOut,
            "overlay_ratio_beforeOverUnfolded_allPtBins.png",
            "Before unfolding / unfolded",
            TString::Format("Before / unfolded ratio overlay, R = %.1f", R).Data()
          );

          DrawRatioOverlaySummary(
            ratioTruthVsUnfoldedHists,
            beforeAfterTruthOut,
            "overlay_ratio_truthOverUnfolded_allPtBins.png",
            "Truth MC / unfolded data",
            TString::Format("Truth / unfolded ratio overlay, R = %.1f", R).Data()
          );

          for (auto* h : perPhoHists) if (h) delete h;
          for (auto* h : perPhoHists_cov) if (h) delete h;
          for (auto* h : perPhoErrRatio) if (h) delete h;
          for (auto* h : perPhoBeforeDataHists) if (h) delete h;
          for (auto* h : perPhoTruthHists) if (h) delete h;
          for (auto* h : perPhoHists_jetEffCorr) if (h) delete h;
          for (auto* h : ratioBeforeVsAfterHists) if (h) delete h;
          for (auto* h : ratioTruthVsUnfoldedHists) if (h) delete h;

          if (hUnfoldTruthGlob) delete hUnfoldTruthGlob;
          if (hUnfoldTruthGlob_cov) delete hUnfoldTruthGlob_cov;
          if (h2JetEff) delete h2JetEff;

          delete hMeasSimGlob;
          delete hTruthSimGlob;
          delete hMeasDataGlob;
          delete hRsp_measXtruth;
          if (h2UnfoldTruth) delete h2UnfoldTruth;
          if (h2UnfoldTruth_cov) delete h2UnfoldTruth_cov;
          if (h2UnfoldTruth_jetEffCorr) delete h2UnfoldTruth_jetEffCorr;

        delete h2RecoData;
        delete h2RecoSim;
        delete h2TruthSim;
        delete h2RspSim;
      }

      // ----------------------------------------------------------------------
      // (B.2) all-radii overlay table (r02, r04, r06) saved to <outBase>/
      //   <triggerName>/unfolding/radii/table2x3_unfolded_perPhoton_dNdXJ_overlay_radii.png
      // ----------------------------------------------------------------------
      {
          auto it02 = perPhoHistsByRKey.find("r02");
          auto it04 = perPhoHistsByRKey.find("r04");
          auto it06 = perPhoHistsByRKey.find("r06");

          const bool have02 = (it02 != perPhoHistsByRKey.end());
          const bool have04 = (it04 != perPhoHistsByRKey.end());
          const bool have06 = (it06 != perPhoHistsByRKey.end());

          bool any = false;
          if (have02) for (auto* h : it02->second) if (h) { any = true; break; }
          if (!any && have04) for (auto* h : it04->second) if (h) { any = true; break; }
          if (!any && have06) for (auto* h : it06->second) if (h) { any = true; break; }

          if (any)
          {
              const auto& analysisRecoBins = UnfoldAnalysisRecoPtBins();
              const int nPtAll = (int)analysisRecoBins.size();

              auto DrawPerPhotonOverlayPad =
                [&](int i)->void
              {
                if (i < 0 || i >= nPtAll)
                {
                  TH1F frame("frame","", 1, 0.0, 2.0);
                  frame.SetMinimum(0.0);
                  frame.SetMaximum(1.0);
                  frame.SetTitle("");
                  frame.GetXaxis()->SetTitle("x_{J}");
                  frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                  frame.Draw("axis");

                  TLatex tx;
                  tx.SetNDC();
                  tx.SetTextFont(42);
                  tx.SetTextSize(0.050);
                  tx.DrawLatex(0.16, 0.50, "MISSING");
                  return;
                }

                const PtBin& b = analysisRecoBins[i];

                TH1* h02 = (have02 && (int)it02->second.size() > i) ? it02->second[i] : nullptr;
                TH1* h04 = (have04 && (int)it04->second.size() > i) ? it04->second[i] : nullptr;
                TH1* h06 = (have06 && (int)it06->second.size() > i) ? it06->second[i] : nullptr;

                double maxY = 0.0;
                auto scanMax = [&](TH1* h)
                {
                    if (!h) return;
                    const int nxb = h->GetNbinsX();
                    for (int ib = 1; ib <= nxb; ++ib)
                    {
                      const double y  = h->GetBinContent(ib);
                      const double ey = h->GetBinError(ib);

                      // Avoid huge "empty-bin" errorbars blowing up the y-range
                      if (y <= 0.0 && ey <= 0.0) continue;

                      const double v = y + ey;
                      if (v > maxY) maxY = v;
                    }
                };
                scanMax(h02);
                scanMax(h04);
                scanMax(h06);

                // Use one of the actual hists to set axes/range (keeps ROOT styling consistent)
                TH1* hBase = (h04 ? h04 : (h02 ? h02 : h06));
                if (hBase)
                {
                    hBase->SetTitle("");
                    hBase->GetXaxis()->SetTitle("x_{J}");
                    hBase->GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                    hBase->SetMinimum(0.0);
                    hBase->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                    hBase->Draw("axis");
                  }
                  else
                  {
                    TH1F frame("frame","", 1, 0.0, 2.0);
                    frame.SetMinimum(0.0);
                    frame.SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                    frame.SetTitle("");
                    frame.GetXaxis()->SetTitle("x_{J}");
                    frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                    frame.Draw("axis");
                }

                if (h02)
                {
                    h02->SetLineColor(kRed);
                    h02->SetMarkerColor(kRed);
                    h02->Draw("E1 X0 same");
                }
                if (h04)
                {
                    h04->SetLineColor(kBlue);
                    h04->SetMarkerColor(kBlue);
                    h04->Draw("E1 X0 same");
                }
                if (h06)
                {
                    h06->SetLineColor(kGreen + 2);
                    h06->SetMarkerColor(kGreen + 2);
                    h06->Draw("E1 X0 same");
                }

                // Centered per-pad title
                {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(22);
                    tx.SetTextSize(0.042);

                    tx.DrawLatex(0.52, 0.955,
                                     TString::Format("Per-photon particle-level x_{J#gamma}, p_{T}^{#gamma} %d-%d GeV",
                                                     b.lo, b.hi).Data());
                }

                // Per-pad cut/trigger label (top-right, 3 lines)
                {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(31);
                    tx.SetTextSize(0.034);

                    const double xR = 0.94;
                    tx.DrawLatex(xR, 0.73, "z_{vtx} < 60 cm");
                    tx.DrawLatex(xR, 0.78, "#Delta #phi > 7#pi/8");
                    tx.DrawLatex(xR, 0.83, "p_{T}^{min, jet} > 5");
                    tx.DrawLatex(xR, 0.88, "Trigger = Photon 4 + MBD NS #geq 1");
                }

                // Legend under the TLatex block (top-right) with VERTICAL-only error bars (EX=0)
                {
                    auto* leg = new TLegend(0.62, 0.5, 0.95, 0.7);
                    leg->SetBorderSize(0);
                    leg->SetFillStyle(0);
                    leg->SetTextFont(42);
                    leg->SetTextSize(0.034);

                    auto makeVertErrLegend = [&](TH1* h) -> TGraphErrors*
                    {
                      if (!h) return nullptr;

                      int ib0 = -1;
                      for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
                      {
                        if (h->GetBinContent(ib) != 0.0 || h->GetBinError(ib) != 0.0) { ib0 = ib; break; }
                      }
                      if (ib0 < 0) ib0 = 1;

                      const double lx  = h->GetXaxis()->GetBinCenter(ib0);
                      const double ly  = h->GetBinContent(ib0);
                      double ley = h->GetBinError(ib0);
                      if (ley <= 0.0) ley = 1.0;

                      auto* g = new TGraphErrors(1);
                      g->SetPoint(0, lx, ly);
                      g->SetPointError(0, 0.0, ley);
                      g->SetMarkerStyle(h->GetMarkerStyle());
                      g->SetMarkerSize(h->GetMarkerSize());
                      g->SetMarkerColor(h->GetMarkerColor());
                      g->SetLineColor(h->GetLineColor());
                      g->SetLineWidth(h->GetLineWidth());
                      return g;
                    };

                    TGraphErrors* gLeg02 = makeVertErrLegend(h02);
                    TGraphErrors* gLeg04 = makeVertErrLegend(h04);
                    TGraphErrors* gLeg06 = makeVertErrLegend(h06);

                    if (gLeg02) leg->AddEntry(gLeg02, "R = 0.2", "pe");
                    if (gLeg04) leg->AddEntry(gLeg04, "R = 0.4", "pe");
                    if (gLeg06) leg->AddEntry(gLeg06, "R = 0.6", "pe");

                    leg->Draw();
                }
              };

              TCanvas c("c_tbl_unf_perPho_overlay_radii", "c_tbl_unf_perPho_overlay_radii", 1800, 1300);
              c.Divide(3, 3, 0.001, 0.001);

              const int nPads = nPtAll;

              for (int ipad = 0; ipad < nPads; ++ipad)
              {
                c.cd(ipad + 1);
                gPad->SetLeftMargin(0.12);
                gPad->SetRightMargin(0.04);
                gPad->SetBottomMargin(0.12);
                gPad->SetTopMargin(0.06);

                DrawPerPhotonOverlayPad(ipad);
              }

              SaveCanvas(c, JoinPath(outBase, "table3x3_unfolded_perPhoton_dNdXJ_overlay_radii.png"));

              for (int i = 0; i < nPtAll; ++i)
              {
                const PtBin& b = analysisRecoBins[i];

                TCanvas cSingle(
                  TString::Format("c_unf_perPho_overlay_radii_%s", b.folder.c_str()).Data(),
                  "c_unf_perPho_overlay_radii_single",
                  900, 700
                );
                ApplyCanvasMargins1D(cSingle);

                DrawPerPhotonOverlayPad(i);

                SaveCanvas(
                  cSingle,
                  JoinPath(
                    outBase,
                    TString::Format("unfolded_perPhoton_dNdXJ_overlay_radii_%s.png", b.folder.c_str()).Data()
                  )
                );
              }
          }
        }

        for (auto& kv : perPhoHistsByRKey)
        {
          for (auto* h : kv.second) if (h) delete h;
        }
        perPhoHistsByRKey.clear();

        delete hPhoRecoData;
        delete hPhoRecoSim;
        delete hPhoTruthSim;
        delete hPhoRespSim;
        delete hPhoResp_measXtruth;
        if (hPhoUnfoldTruth) delete hPhoUnfoldTruth;
        if (hPhoUnfoldTruth_cov) delete hPhoUnfoldTruth_cov;

        if (gAtlasPP) delete gAtlasPP;
    }

    void RunPurityCorrectedUncorrectedOverlayPP(Dataset& dsData, Dataset& dsSim)
    {
      const string simDataBase = JoinPath(dsSim.outBase, dsData.trigger);
      const string inBaseUnc = JoinPath(simDataBase, "unfolding/nonPurityCorrected/radii");
      const string inBaseCor = JoinPath(simDataBase, "unfolding/purityCorrected/radii");
      const string outBase   = JoinPath(simDataBase, "unfolding/purityCorrectedUncorrectedOverly");

      EnsureDir(outBase);

      auto FindStoredUnfoldPtIndexForAnalysis = [&](const PtBin& b)->int
        {
          const auto& analysisBins = UnfoldAnalysisRecoPtBins();
          for (int i = 0; i < (int)analysisBins.size(); ++i)
          {
            if (analysisBins[i].lo == b.lo && analysisBins[i].hi == b.hi) return i;
          }
          return -1;
      };

      for (const auto& rKey : kRKeys)
      {
        const string uncRoot = JoinPath(JoinPath(inBaseUnc, rKey), "rooUnfold_outputs.root");
        const string corRoot = JoinPath(JoinPath(inBaseCor, rKey), "rooUnfold_outputs.root");
        const string rOut    = JoinPath(outBase, rKey);

        EnsureDir(rOut);

        TFile* fUnc = TFile::Open(uncRoot.c_str(), "READ");
        TFile* fCor = TFile::Open(corRoot.c_str(), "READ");

        if (!fUnc || !fCor || fUnc->IsZombie() || fCor->IsZombie())
        {
          if (fUnc) { fUnc->Close(); delete fUnc; }
          if (fCor) { fCor->Close(); delete fCor; }
          continue;
        }

        auto GetH = [&](TFile* f, int iStored, const string& stem)->TH1*
        {
            if (!f)
            {
              cout << ANSI_BOLD_RED
                   << "[PURITY OVERLAY DEBUG] null TFile for stem=" << stem
                   << "  iStored=" << iStored
                   << ANSI_RESET << "\n";
              return nullptr;
            }

            if (iStored < 0)
            {
              cout << ANSI_BOLD_YEL
                   << "[PURITY OVERLAY DEBUG] invalid stored index for stem=" << stem
                   << "  iStored=" << iStored
                   << ANSI_RESET << "\n";
              return nullptr;
            }

            const string hname = TString::Format("%s_pTbin%d", stem.c_str(), iStored + 1).Data();

            cout << ANSI_BOLD_CYN
                 << "[PURITY OVERLAY DEBUG] file=" << f->GetName()
                 << "  requesting hist=" << hname
                 << ANSI_RESET << "\n";

            TH1* h = dynamic_cast<TH1*>(f->Get(hname.c_str()));
            if (!h)
            {
              cout << ANSI_BOLD_RED
                   << "  -> missing histogram: " << hname
                   << ANSI_RESET << "\n";
              return nullptr;
            }

            TH1* hc = dynamic_cast<TH1*>(h->Clone(TString::Format("%s_clone_%d", hname.c_str(), iStored + 1).Data()));
            if (!hc)
            {
              cout << ANSI_BOLD_RED
                   << "  -> clone failed for histogram: " << hname
                   << ANSI_RESET << "\n";
              return nullptr;
            }

            hc->SetDirectory(nullptr);
            EnsureSumw2(hc);

            cout << "  -> found " << hc->GetName()
                 << "  nbins=" << hc->GetNbinsX()
                 << "  integral=" << hc->Integral(1, hc->GetNbinsX())
                 << "  integral(width)=" << hc->Integral(1, hc->GetNbinsX(), "width")
                 << "\n";

            for (int ib = 0; ib <= hc->GetNbinsX() + 1; ++ib)
            {
              if (ib == 0)
              {
                cout << "     UF  content=" << hc->GetBinContent(ib)
                     << "  error=" << hc->GetBinError(ib) << "\n";
              }
              else if (ib == hc->GetNbinsX() + 1)
              {
                cout << "     OF  content=" << hc->GetBinContent(ib)
                     << "  error=" << hc->GetBinError(ib) << "\n";
              }
              else
              {
                cout << "     bin " << ib
                     << "  [" << hc->GetXaxis()->GetBinLowEdge(ib)
                     << ","   << hc->GetXaxis()->GetBinUpEdge(ib) << "]"
                     << "  content=" << hc->GetBinContent(ib)
                     << "  error="   << hc->GetBinError(ib) << "\n";
              }
            }

            return hc;
        };

          auto DrawPad = [&](int iAnalysis)->void
          {
              const auto& analysisBins = UnfoldAnalysisRecoPtBins();

              if (iAnalysis < 0 || iAnalysis >= (int)analysisBins.size())
              {
                TH1F frame("frame","", 1, 0.0, 2.0);
                frame.SetMinimum(0.0);
                frame.SetMaximum(1.0);
                frame.SetTitle("");
                frame.GetXaxis()->SetTitle("x_{J}");
                frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                frame.Draw("axis");
                return;
            }

            const PtBin& b = analysisBins[iAnalysis];
            const int iStored = FindStoredUnfoldPtIndexForAnalysis(b);

            TH1* hUnc = GetH(fUnc, iStored, "h_xJ_unf_perPho");
            TH1* hCor = GetH(fCor, iStored, "h_xJ_unf_perPho");

            auto HasFiniteDrawableContent = [&](TH1* h)->bool
            {
              if (!h) return false;
              const int nb = h->GetNbinsX();
              for (int ib = 1; ib <= nb; ++ib)
              {
                const double y  = h->GetBinContent(ib);
                const double ey = h->GetBinError(ib);
                if (std::isfinite(y) && std::isfinite(ey) && (y != 0.0 || ey != 0.0)) return true;
              }
              return false;
            };

            auto SanitizeHistForDraw = [&](TH1* h)->void
            {
              if (!h) return;
              const int nb = h->GetNbinsX();
              for (int ib = 0; ib <= nb + 1; ++ib)
              {
                double y  = h->GetBinContent(ib);
                double ey = h->GetBinError(ib);
                if (!std::isfinite(y))  y  = 0.0;
                if (!std::isfinite(ey)) ey = 0.0;
                h->SetBinContent(ib, y);
                h->SetBinError  (ib, ey);
              }
            };

            SanitizeHistForDraw(hUnc);
            SanitizeHistForDraw(hCor);

            const bool haveUncFinite = HasFiniteDrawableContent(hUnc);
            const bool haveCorFinite = HasFiniteDrawableContent(hCor);

            double maxY = 0.0;
            auto scan = [&](TH1* h)
            {
              if (!h) return;
              h->GetXaxis()->SetRangeUser(0.0, 2.0);
              const int nb = h->GetNbinsX();
              for (int ib = 1; ib <= nb; ++ib)
              {
                const double y  = h->GetBinContent(ib);
                const double ey = h->GetBinError(ib);
                if (!std::isfinite(y) || !std::isfinite(ey)) continue;
                if (y <= 0.0 && ey <= 0.0) continue;
                maxY = std::max(maxY, y + ey);
              }
            };

            if (haveUncFinite) scan(hUnc);
            if (haveCorFinite) scan(hCor);

            TH1* hBase = (haveUncFinite ? hUnc : (haveCorFinite ? hCor : nullptr));
            TH1* hUncDraw = nullptr;
            TH1* hCorDraw = nullptr;

            if (hBase)
            {
                hBase->GetXaxis()->SetRangeUser(0.0, 2.0);
                hBase->SetTitle("");
                hBase->GetXaxis()->SetTitle("x_{J}");
                hBase->GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                hBase->SetMinimum(0.0);
                hBase->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                hBase->DrawClone("axis");
            }
            else
              {
                TH1F frame("frame","", 1, 0.0, 2.0);
                frame.SetMinimum(0.0);
                frame.SetMaximum(1.0);
                frame.SetTitle("");
                frame.GetXaxis()->SetTitle("x_{J}");
                frame.GetYaxis()->SetTitle("(1/N_{#gamma}) dN/dx_{J}");
                frame.DrawCopy("axis");
            }

            if (haveUncFinite)
            {
                hUnc->SetLineColor(kBlack);
                hUnc->SetMarkerColor(kBlack);
                hUnc->SetMarkerStyle(20);
                hUnc->SetMarkerSize(1.0);
                hUnc->SetLineWidth(2);
                hUncDraw = dynamic_cast<TH1*>(hUnc->DrawClone("E1 X0 same"));
            }

            if (haveCorFinite)
            {
                hCor->SetLineColor(kBlue + 1);
                hCor->SetMarkerColor(kBlue + 1);
                hCor->SetMarkerStyle(24);
                hCor->SetMarkerSize(1.0);
                hCor->SetLineWidth(2);
                hCorDraw = dynamic_cast<TH1*>(hCor->DrawClone("E1 X0 same"));
            }

            TLatex tx;
            tx.SetNDC();
            tx.SetTextFont(42);
            tx.SetTextAlign(22);
            tx.SetTextSize(0.042);
            tx.DrawLatex(0.50, 0.955,
                           TString::Format("Per-photon particle-level x_{J#gamma}, p_{T}^{#gamma} %d-%d GeV",
                                           b.lo, b.hi).Data());

              if (haveUncFinite || haveCorFinite)
              {
                  auto MakeVerticalErrLegendProxy = [&](TH1* hSrc, const char* gname) -> TGraphErrors*
                  {
                    if (!hSrc) return nullptr;

                    int ib0 = -1;
                    for (int ib = 1; ib <= hSrc->GetNbinsX(); ++ib)
                    {
                      const double y  = hSrc->GetBinContent(ib);
                      const double ey = hSrc->GetBinError(ib);
                      if (std::isfinite(y) && std::isfinite(ey) && (y != 0.0 || ey != 0.0))
                      {
                        ib0 = ib;
                        break;
                      }
                    }
                    if (ib0 < 0) ib0 = 1;

                    const double x  = hSrc->GetXaxis()->GetBinCenter(ib0);
                    const double y  = hSrc->GetBinContent(ib0);
                    double ey = hSrc->GetBinError(ib0);
                    if (!std::isfinite(ey) || ey <= 0.0) ey = 1.0;

                    TGraphErrors* g = new TGraphErrors(1);
                    g->SetName(gname);
                    g->SetPoint(0, x, y);
                    g->SetPointError(0, 0.0, ey);
                    g->SetMarkerStyle(hSrc->GetMarkerStyle());
                    g->SetMarkerSize(hSrc->GetMarkerSize());
                    g->SetMarkerColor(hSrc->GetMarkerColor());
                    g->SetLineColor(hSrc->GetLineColor());
                    g->SetLineWidth(hSrc->GetLineWidth());
                    return g;
                  };

                  TGraphErrors* gLegUnc = MakeVerticalErrLegendProxy(hUncDraw, "gLegUnc_purityOverlay");
                  TGraphErrors* gLegCor = MakeVerticalErrLegendProxy(hCorDraw, "gLegCor_purityOverlay");

                  TLegend leg(0.58, 0.72, 0.94, 0.88);
                  leg.SetBorderSize(0);
                  leg.SetFillStyle(0);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.032);
                  if (gLegUnc) leg.AddEntry(gLegUnc, "non-purity-corrected", "pe");
                  if (gLegCor) leg.AddEntry(gLegCor, "purity-corrected", "pe");
                  leg.DrawClone();

                  if (gLegUnc) delete gLegUnc;
                  if (gLegCor) delete gLegCor;
              }

            if (gPad) { gPad->Modified(); gPad->Update(); }

            if (hUnc) delete hUnc;
            if (hCor) delete hCor;
          };

        TCanvas c(TString::Format("c_tbl_purityOv_%s", rKey.c_str()).Data(), "c_tbl_purityOv", 1800, 1300);
        c.Divide(3, 3, 0.001, 0.001);

        const auto& analysisBins = UnfoldAnalysisRecoPtBins();
        const int nPads = (int)analysisBins.size();

        for (int ipad = 0; ipad < nPads; ++ipad)
        {
              c.cd(ipad + 1);
              gPad->SetLeftMargin(0.12);
              gPad->SetRightMargin(0.04);
              gPad->SetBottomMargin(0.12);
              gPad->SetTopMargin(0.06);

              DrawPad(ipad);
        }

        SaveCanvas(c, JoinPath(rOut, "table3x3_purityCorrected_vs_nonPurityCorrected_perPhoton_dNdXJ.png"));

        fUnc->Close();
        fCor->Close();
        delete fUnc;
        delete fCor;
      }
    }
    // =============================================================================
    // Au+Au SIM+DATA unfolding wrapper
    //
    // Matches per-centrality SIM and DATA datasets, then for each pair runs
    // the standard unfold pipeline with all 4 (purity × combinatoric) variants.
    //
    // Output structure (under each centrality SIM outBase):
    //   <trigAA>/unfolding/{purityCorrected,nonPurityCorrected}/
    //       {combinatoricSub,nonCombSub}/radii/<rKey>/...
    // =============================================================================
    void RunRooUnfoldPipeline_SimAndDataAUAU(vector<Dataset>& datasets)
    {
      cout << ANSI_BOLD_CYN << "\n==============================\n"
           << "[SECTION 5I-AA] RooUnfold pipeline (SIM+DATA AuAu)\n"
           << "==============================" << ANSI_RESET << "\n";

      const bool oldPurity = gApplyPurityCorrectionForUnfolding;
      const bool oldComb   = gApplyCombinatoricSubtractionForUnfolding;

      // --- Match each DATA trigger+centrality dataset to the SIM centrality dataset ---
      map<string, Dataset*> simByCent;
      vector<Dataset*> dataDatasets;

      for (auto& ds : datasets)
      {
        if (ds.isSim) simByCent[ds.centSuffix] = &ds;
        else          dataDatasets.push_back(&ds);
      }

      struct Variant
      {
        bool purity;
        bool comb;
        const char* folderLabel;
      };

      const vector<Variant> variants = {
        {false, false, "nonPurityCorrected"},
        {true,  false, "purityCorrected"},
        {false, true,  "combinatoricJetSubtractectedNoPurityCorrection"},
        {true,  true,  "combinatoricJetSubtractectedWithPurityCorrection"}
      };

      int nPairs = 0;
      for (auto* dsDATA : dataDatasets)
      {
        auto it = simByCent.find(dsDATA->centSuffix);
        if (it == simByCent.end())
        {
          cout << ANSI_BOLD_YEL
               << "[WARN] No SIM dataset matches DATA trigger '" << dsDATA->trigger
               << "' centrality '" << dsDATA->centSuffix << "'. Skipping.\n"
               << ANSI_RESET;
          continue;
        }
        Dataset* dsSIM = it->second;

        cout << ANSI_BOLD_CYN
             << "\n--- AuAu unfold pair: " << dsSIM->centLabel
             << " | trigger=" << dsDATA->trigger << " ---\n"
             << "  SIM:  " << dsSIM->label  << "  outBase=" << dsSIM->outBase << "\n"
             << "  DATA: " << dsDATA->label << "  outBase=" << dsDATA->outBase << "\n"
             << ANSI_RESET;

        for (const auto& v : variants)
        {
          gApplyPurityCorrectionForUnfolding        = v.purity;
          gApplyCombinatoricSubtractionForUnfolding = v.comb;

          cout << "  -> [AuAu unfold] " << dsSIM->centLabel
               << " | trigger=" << dsDATA->trigger
               << " | outputFolder=" << v.folderLabel << " ...\n";

          RunRooUnfoldPipeline_SimAndDataPP(*dsDATA, *dsSIM);

          cout << "     [OK] " << v.folderLabel << "\n";
        }

        ++nPairs;
      }

      gApplyPurityCorrectionForUnfolding        = oldPurity;
      gApplyCombinatoricSubtractionForUnfolding = oldComb;

      cout << ANSI_BOLD_CYN
           << "[DONE] AuAu unfold pipeline: processed " << nPairs << " trigger-centrality pairs.\n"
           << ANSI_RESET;
    }
#else
    void RunRooUnfoldPipeline_SimAndDataPP(Dataset& dsData, Dataset& dsSim)
    {
      (void)dsData;
      (void)dsSim;
      cout << ANSI_BOLD_YEL
           << "[WARN] RooUnfold headers not found at compile time (ARJ_HAVE_ROOUNFOLD=0). Skipping RooUnfold pipeline."
           << ANSI_RESET << "\n";
    }

    void RunRooUnfoldPipeline_SimAndDataAUAU(vector<Dataset>& datasets)
    {
      (void)datasets;
      cout << ANSI_BOLD_YEL
           << "[WARN] RooUnfold headers not found at compile time (ARJ_HAVE_ROOUNFOLD=0). Skipping AuAu RooUnfold pipeline."
           << ANSI_RESET << "\n";
    }

    void RunPurityCorrectedUncorrectedOverlayPP(Dataset& dsData, Dataset& dsSim)
    {
      (void)dsData;
      (void)dsSim;
      cout << ANSI_BOLD_YEL
           << "[WARN] RooUnfold headers not found at compile time (ARJ_HAVE_ROOUNFOLD=0). Skipping purity-corrected vs non-purity-corrected overlay."
           << ANSI_RESET << "\n";
    }
#endif
