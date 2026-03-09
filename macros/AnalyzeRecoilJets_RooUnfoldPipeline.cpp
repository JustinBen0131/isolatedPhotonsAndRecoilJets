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
        const string outBase = JoinPath(dsData.outBase, "unfolding/" + unfoldVariant + "/radii");
        int mkrc = -999;
        if (gSystem) mkrc = gSystem->mkdir(outBase.c_str(), true);
        else EnsureDir(outBase);

        cout << "  Unfolding variant: " << unfoldVariant << "\n";
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

        const int kBayesIterPho = 3;

        // Toy settings:
        //   - "final" is for the baseline unfolded spectrum that feeds your main outputs
        //   - "scan" is for iteration-stability / closure / half-closure diagnostics
        const int kNToysPhoFinal = 600;
        const int kNToysPhoScan  = 120;

        RooUnfoldResponse respPho(hPhoRecoSim, hPhoTruthSim, hPhoResp_measXtruth, "respPho", "respPho");

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

                  TLegend leg(0.55, 0.74, 0.88, 0.88);
                  leg.SetBorderSize(0);
                  leg.SetFillStyle(0);
                  leg.SetTextFont(42);
                  leg.SetTextSize(0.035);
                  leg.AddEntry(&gLegData, "#epsilon_{#gamma}^{eff,data} (reco / unfolded truth)", "pe");
                  leg.AddEntry(&gLegSim,  "#epsilon_{#gamma}^{MC} (1 - misses/truth)", "pe");
                  leg.Draw();

                  drawLineAtOne(hEffSim);

                  DrawLatexLines(0.14,0.92, DefaultHeaderLines(dsData), 0.034, 0.045);
                  DrawLatexLines(0.14,0.78, { "Photon efficiency diagnostics: DATA implied vs SIM truth-matching" }, 0.030, 0.040);

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
          const int kMaxIt = 10;

          vector<double> xIt, exIt, yRelStat, eyRelStat, yRelDev, eyRelDev;

          TH1* hPrev = nullptr;

            for (int it = 1; it <= kMaxIt; ++it)
            {
              RooUnfoldBayes uIt(&respPho, hPhoRecoData, it);
              uIt.SetVerbose(0);
              uIt.SetNToys(kNToysPhoScan);

              TH1* hIt = nullptr;
              if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
              hIt = uIt.Hreco(RooUnfold::kCovToy);
              if (gSystem) gSystem->RedirectOutput(0);
              if (!hIt) continue;

              hIt->SetDirectory(nullptr);
              EnsureSumw2(hIt);

              const int nb = hIt->GetNbinsX();

              double sumV  = 0.0;
              double sumE2 = 0.0;

              for (int ib = 1; ib <= nb; ++ib)
              {
                const double v  = hIt->GetBinContent(ib);
                const double ev = hIt->GetBinError(ib);
                if (v == 0.0 && ev == 0.0) continue;
                sumV  += v;
                sumE2 += ev * ev;
              }

              double relStat = 0.0;
              if (sumV > 0.0) relStat = std::sqrt(sumE2) / sumV;

              if (!hPrev)
              {
                hPrev = hIt;
                continue;
              }

              double relDev = 0.0;
              {
                double num2 = 0.0;
                double den2 = 0.0;

                for (int ib = 1; ib <= nb; ++ib)
                {
                  const double v  = hIt->GetBinContent(ib);
                  const double vp = hPrev->GetBinContent(ib);
                  const double d  = v - vp;

                  num2 += d * d;
                  den2 += v * v;
                }

                if (den2 > 0.0) relDev = std::sqrt(num2 / den2);
              }

              xIt.push_back((double)it);
              exIt.push_back(0.0);
              yRelStat.push_back(relStat);
              eyRelStat.push_back(0.0);
              yRelDev.push_back(relDev);
              eyRelDev.push_back(0.0);

              delete hPrev;
              hPrev = hIt;
            }

          if (hPrev) delete hPrev;

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

            // total relative deviation (it vs it-1) (BLUE)
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
                hRat->GetYaxis()->SetTitle("Closure: unfolded MC / truth MC");
                hRat->SetMarkerStyle(24);
                hRat->SetMarkerSize(1.1);
                hRat->SetLineWidth(2);

                TCanvas c("c_pho_closure", "c_pho_closure", 900, 700);
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

                  // Top-left Bayes + Step label (match 2D closure styling)
                  {
                    TLatex tx;
                    tx.SetNDC();
                    tx.SetTextFont(42);
                    tx.SetTextAlign(13);
                    tx.SetTextSize(0.038);
                    tx.DrawLatex(0.15, 0.90, TString::Format("Bayes it=%d (photon)", kBayesIterPho).Data());
                    tx.DrawLatex(0.15, 0.855, "1D photon-yield unfolding");
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
                      tx.DrawLatex(0.15, 0.90, TString::Format("Bayes it=%d (photon)", kBayesIterPho).Data());
                      tx.DrawLatex(0.15, 0.855, "1D photon-yield unfolding");
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

      // -------------------------------------------------------------------------
      // (B) Jet unfolding per radius: unfold global-bin vector and unflatten back
      // -------------------------------------------------------------------------
      const int kBayesIterXJ = 3;

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
        EnsureDir(beforeAfterDataOut);
        EnsureDir(beforeAfterTruthOut);

          const string nameReco   = "h2_unfoldReco_pTgamma_xJ_incl_"           + rKey;
          const string nameRecoC  = "h2_unfoldReco_pTgamma_xJ_incl_sidebandC_" + rKey;
          const string nameTruth  = "h2_unfoldTruth_pTgamma_xJ_incl_"          + rKey;
          const string nameRsp    = "h2_unfoldResponse_pTgamma_xJ_incl_"       + rKey;

          TH2* h2RecoData_in       = GetObj<TH2>(dsData, nameReco,  true, true, true);
          TH2* h2RecoData_sideC_in = (gApplyPurityCorrectionForUnfolding
                                      ? GetObj<TH2>(dsData, nameRecoC, true, true, true)
                                      : nullptr);
          TH2* h2RecoSim_in        = GetObj<TH2>(dsSim,  nameReco,  true, true, true);
          TH2* h2TruthSim_in       = GetObj<TH2>(dsSim,  nameTruth, true, true, true);
          TH2* h2RspSim_in         = GetObj<TH2>(dsSim,  nameRsp,   true, true, true);

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

            // required for log-z (must be > 0)
            hScat->SetMinimum(0.5);

            hScat->Draw("COLZ");

            {
                TLatex tx;
                tx.SetNDC();
                tx.SetTextFont(42);
                tx.SetTextAlign(22);
                tx.SetTextSize(0.040);
                tx.DrawLatex(0.50, 0.965,
                             TString::Format("2D Response Matrix, Photon 10 + 20 GeV merged Pythia, R = %.1f", R).Data());
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
        vector<TH1*> ratioBeforeVsAfterHists(nPtAll, nullptr);
        vector<TH1*> ratioTruthVsUnfoldedHists(nPtAll, nullptr);

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
            TH1* hMeasSimGlob_closure = CloneTH1(hMeasSimGlob,
              TString::Format("hMeasSimGlob_forClosure_%s", rKey.c_str()).Data());
            if (hMeasSimGlob_closure) hMeasSimGlob_closure->SetDirectory(nullptr);

            RooUnfoldBayes unfoldXJ_closure(&respXJ, (hMeasSimGlob_closure ? hMeasSimGlob_closure : hMeasSimGlob), kBayesIterXJ);
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

            vector<double> xPt, exPt, yRat, eyRat;

            const auto& analysisRecoBins = UnfoldAnalysisRecoPtBins();
            const int nPtClosure = (int)analysisRecoBins.size();
            for (int i = 0; i < nPtClosure; ++i)
            {
              const PtBin& b = analysisRecoBins[i];
              const double cen = 0.5 * (b.lo + b.hi);
              const double ex  = 0.5 * (b.hi - b.lo);

              if (!h2UnfoldTruth_closure || !h2TruthSim) continue;

              const int ixTruth = h2TruthSim->GetXaxis()->FindBin(cen);
              if (ixTruth < 1 || ixTruth > h2TruthSim->GetXaxis()->GetNbins()) continue;

              TH1D* hXJ_truth = h2TruthSim->ProjectionY(
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
                  tx.DrawLatex(0.15, 0.90, TString::Format("Bayes it=%d (xJ)", kBayesIterXJ).Data());
                  tx.DrawLatex(0.15, 0.855, "2D (p_{T}^{#gamma}, x_{J}) unfolding");
                }

                SaveCanvas(c, JoinPath(rOut, "closure_unfoldedOverTruth_integral_vs_pTgamma.png"));
            }

            if (h2UnfoldTruth_closure) delete h2UnfoldTruth_closure;
            if (hUnfoldTruthGlob_closure) delete hUnfoldTruthGlob_closure;
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
                      tx.DrawLatex(0.15, 0.90, TString::Format("Bayes it=%d (xJ)", kBayesIterXJ).Data());
                      tx.DrawLatex(0.15, 0.855, "2D (p_{T}^{#gamma}, x_{J}) unfolding");
                    }

                    SaveCanvas(c, JoinPath(rOut, "halfClosure_unfoldedOverTruth_integral_vs_pTgamma.png"));
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

                SaveCanvas(c, JoinPath(rOut, TString::Format("xJ_unfolded_perPhoton_pTbin%d.png", i + 1).Data()));

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
                // -------------------------------------------------------------------
                if (gAtlasPP)
                {
                  TCanvas cO(TString::Format("c_perPho_LHC_%s_%d", rKey.c_str(), i + 1).Data(), "c_perPho_LHC", 900, 700);
                  ApplyCanvasMargins1D(cO);

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
                    hTmp->SetMaximum((maxY > 0.0) ? (1.15 * maxY) : 1.0);
                    hTmp->GetXaxis()->SetRangeUser(0.0, 2.0);

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

                      // Legend: top-left
                      TLegend leg(0.53,0.76,0.85,0.90);
                      leg.SetTextFont(42);
                      leg.SetTextSize(0.029);
                      leg.AddEntry(&gSphLeg,
                                     TString::Format("sPHENIX unfolded, p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data(),
                                     "pe");
                      leg.AddEntry(gAtlasPP,
                                     TString::Format("ATLAS unfolded, p_{T}^{#gamma} = %s", kAtlasTable1PhoPtLabel.c_str()).Data(),
                                     "pe");
                      leg.Draw();

                    SaveCanvas(cO, JoinPath(overlayOut, TString::Format("xJ_unfolded_perPhoton_LHCoverlay_pTbin%d.png", i + 1).Data()));

                    delete hTmp;
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

              SaveCanvas(c, JoinPath(rOut, "table3x3_unfolded_perPhoton_dNdXJ.png"));

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
          // Output: <rOut>/unfold_iterStability_relChange_relStat.png
          // ----------------------------------------------------------------------
          {
            const int kMaxIterPlot = 10;

            cout << ANSI_BOLD_CYN
                 << "[UNF ITER QA] Building iteration-stability plot\n"
                 << "  rKey=" << rKey << "  R=" << std::fixed << std::setprecision(1) << R << "\n"
                 << "  outdir=" << rOut << "\n"
                 << "  kMaxIterPlot=" << kMaxIterPlot << "\n"
                 << "  respXJ ptr=" << (void*)&respXJ << "  hMeasDataGlob ptr=" << (void*)hMeasDataGlob << "\n"
                 << ANSI_RESET;

            std::vector<double> xIt, exIt, yRelStat, eyRelStat, yRelChange, eyRelChange;

            TH1* hPrev = nullptr;

            for (int it = 1; it <= kMaxIterPlot; ++it)
            {
              if (!hMeasDataGlob)
              {
                cout << ANSI_BOLD_RED << "[UNF ITER QA][FATAL] hMeasDataGlob is null. Skipping." << ANSI_RESET << "\n";
                break;
              }

              RooUnfoldBayes u(&respXJ, hMeasDataGlob, it);
              u.SetVerbose(0);
              u.SetNToys(kNToysXJScan);

              TH1* hCurr = nullptr;
              if (gSystem) gSystem->RedirectOutput("/dev/null", "w");
              hCurr = u.Hreco(RooUnfold::kCovToy);
              if (gSystem) gSystem->RedirectOutput(0);
              if (!hCurr)
              {
                  cout << ANSI_BOLD_RED << "[UNF ITER QA][WARN] Hreco returned null at it=" << it << ANSI_RESET << "\n";
                  continue;
              }

              hCurr->SetDirectory(nullptr);
              EnsureSumw2(hCurr);

              const int nb = hCurr->GetNbinsX();

              double sumV2 = 0.0;
              double sumE2 = 0.0;

              for (int ib = 1; ib <= nb; ++ib)
              {
                const double v  = hCurr->GetBinContent(ib);
                const double ev = hCurr->GetBinError(ib);
                sumV2 += v * v;
                sumE2 += ev * ev;
              }

              const double relStat = (sumV2 > 0.0) ? std::sqrt(sumE2 / sumV2) : 0.0;

              double relChg = 0.0;
              if (hPrev)
              {
                double sumD2 = 0.0;
                for (int ib = 1; ib <= nb; ++ib)
                {
                  const double d = hCurr->GetBinContent(ib) - hPrev->GetBinContent(ib);
                  sumD2 += d * d;
                }
                relChg = (sumV2 > 0.0) ? std::sqrt(sumD2 / sumV2) : 0.0;
              }

              cout << ANSI_BOLD_YEL
                   << "[UNF ITER QA] it=" << it
                   << "  nb=" << nb
                   << "  relStat=" << std::setprecision(6) << relStat
                   << "  relChange(it vs it-1)=" << std::setprecision(6) << relChg
                   << ANSI_RESET << "\n";

              // Plot starts at it=2 (but uses it=1 internally to compute change at it=2)
              if (it >= 2)
              {
                xIt.push_back((double)it);
                exIt.push_back(0.0);

                yRelStat.push_back(relStat);
                eyRelStat.push_back(0.0);

                yRelChange.push_back(relChg);
                eyRelChange.push_back(0.0);
              }

              if (hPrev) delete hPrev;
              hPrev = hCurr;
            }

            if (hPrev) delete hPrev;

            cout << ANSI_BOLD_CYN
                 << "[UNF ITER QA] Points prepared: n=" << xIt.size()
                 << " (expected " << std::max(0, kMaxIterPlot - 1) << " for it=2..kMax)\n"
                 << ANSI_RESET;

            if (!xIt.empty())
            {
              double yMax = 0.0;
              for (size_t i = 0; i < yRelStat.size();   ++i) yMax = std::max(yMax, yRelStat[i]);
              for (size_t i = 0; i < yRelChange.size(); ++i) yMax = std::max(yMax, yRelChange[i]);
              if (yMax <= 0.0) yMax = 1.0;

              cout << ANSI_BOLD_CYN << "[UNF ITER QA] yMax=" << std::setprecision(6) << yMax << ANSI_RESET << "\n";

              TCanvas cIt(TString::Format("c_iterStability_%s", rKey.c_str()).Data(), "c_iterStability", 900, 700);
              ApplyCanvasMargins1D(cIt);

              // x-axis starts at 1 (but first point is it=2)
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

              // Legend higher in the top-left
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
                  tx.DrawLatex(0.50, 0.965, "Iteration Stability, R = 0.4, Photon 4 + MBD NS #geq 1, Run24pp");
              }

              cIt.Modified();
              cIt.Update();

              const std::string outPng = JoinPath(rOut, "unfold_iterStability_relChange_relStat.png");
              cout << ANSI_BOLD_CYN << "[UNF ITER QA] Saving: " << outPng << ANSI_RESET << "\n";
              SaveCanvas(cIt, outPng);
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
              if (h2UnfoldTruth)        h2UnfoldTruth->Write("h2_truthUnfold_pTgamma_xJ");
              if (h2UnfoldTruth_cov)    h2UnfoldTruth_cov->Write("h2_truthUnfold_pTgamma_xJ_covariance");

              for (int i = 0; i < nPtAll; ++i)
              {
                  if (perPhoHists[i])     perPhoHists[i]->Write(TString::Format("h_xJ_unf_perPho_pTbin%d", i + 1).Data());
                  if (perPhoHists_cov[i]) perPhoHists_cov[i]->Write(TString::Format("h_xJ_unf_perPho_cov_pTbin%d", i + 1).Data());
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
          for (auto* h : ratioBeforeVsAfterHists) if (h) delete h;
          for (auto* h : ratioTruthVsUnfoldedHists) if (h) delete h;

          if (hUnfoldTruthGlob) delete hUnfoldTruthGlob;
          if (hUnfoldTruthGlob_cov) delete hUnfoldTruthGlob_cov;

          delete hMeasSimGlob;
          delete hTruthSimGlob;
          delete hMeasDataGlob;
          delete hRsp_measXtruth;
          if (h2UnfoldTruth) delete h2UnfoldTruth;
          if (h2UnfoldTruth_cov) delete h2UnfoldTruth_cov;

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

    void RunPurityCorrectedUncorrectedOverlayPP(Dataset& dsData)
    {
      const string inBaseUnc = JoinPath(dsData.outBase, "unfolding/nonPurityCorrected/radii");
      const string inBaseCor = JoinPath(dsData.outBase, "unfolding/purityCorrected/radii");
      const string outBase   = JoinPath(dsData.outBase, "unfolding/purityCorrectedUncorrectedOverly");

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
#else
    void RunRooUnfoldPipeline_SimAndDataPP(Dataset& dsData, Dataset& dsSim)
    {
      (void)dsData;
      (void)dsSim;
      cout << ANSI_BOLD_YEL
           << "[WARN] RooUnfold headers not found at compile time (ARJ_HAVE_ROOUNFOLD=0). Skipping RooUnfold pipeline."
           << ANSI_RESET << "\n";
    }

    void RunPurityCorrectedUncorrectedOverlayPP(Dataset& dsData)
    {
      (void)dsData;
      cout << ANSI_BOLD_YEL
           << "[WARN] RooUnfold headers not found at compile time (ARJ_HAVE_ROOUNFOLD=0). Skipping purity-corrected vs non-purity-corrected overlay."
           << ANSI_RESET << "\n";
    }
#endif
