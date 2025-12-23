//======================================================================
//  Fun4All_recoilJets.C
//  --------------------------------------------------------------------
#pragma once
#if defined(__CINT__) || defined(__CLING__)
  R__ADD_INCLUDE_PATH(/sphenix/u/patsfan753/thesisAnalysis/install/include)
#endif
#if defined(__CLING__)
  #pragma cling add_include_path("/sphenix/u/patsfan753/thesisAnalysis/install/include")
#endif
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)

//–––– Standard Fun4All / sPHENIX ––––––––––––––––––––––––––––––––––––––
#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllUtils.h>
#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <calobase/TowerInfoContainer.h>
#include <mbd/MbdPmtContainer.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/CaloWaveformProcessing.h>
#include <caloreco/CaloTowerBuilder.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <calotrigger/MinimumBiasClassifier.h>
#include <ffamodules/FlagHandler.h>
#include <ffamodules/CDBInterface.h>
#include <clusteriso/ClusterIso.h>
#include <calotrigger/TriggerRunInfoReco.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <caloreco/CaloGeomMapping.h>
#include <caloreco/RawClusterPositionCorrection.h>
#include <calobase/RawTowerGeom.h>
#include <caloreco/RawTowerCalibration.h>
#include <caloreco/PhotonClusterBuilder.h>
#include <jetbase/Jet.h>

#include <jetbase/FastJetOptions.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <globalvertex/GlobalVertexReco.h>
#include <caloreco/CaloTowerCalib.h>

#include <centrality/CentralityReco.h>
#include <mbd/MbdEvent.h>
#include <mbd/MbdReco.h>
#include <zdcinfo/ZdcReco.h>
#include <phool/recoConsts.h>
#include <phool/PHRandomSeed.h>
#include <jetbase/FastJetOptions.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include "/sphenix/u/patsfan753/scratch/thesisAnalysis/src/RecoilJets.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <cstdlib>     // getenv
#include <algorithm>   // std::transform
#include <cctype>      // std::tolower
#include <dlfcn.h>     // dlopen RTLD_NOLOAD
#include <TSystem.h>   // gSystem, GetBuildArch/Compiler info
#include <Calo_Calib.C>

// Load your local overrides FIRST so their symbols/dictionaries win
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_reco.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libcalo_io.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libclusteriso.so)
R__LOAD_LIBRARY(/sphenix/u/patsfan753/thesisAnalysis/install/lib/libRecoilJets.so)

// Then load the rest of the environment stack
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffarawobjects.so)
R__LOAD_LIBRARY(libcaloTreeGen.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libglobalvertex.so)
R__LOAD_LIBRARY(libcentrality.so)      // always
R__LOAD_LIBRARY(libcentrality_io.so)   // if you instantiate CentralityReco
R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libzdcinfo.so)
R__LOAD_LIBRARY(libmbd.so)


//======================================================================
//  Convenience helpers
//======================================================================
namespace detail
{
  /// Throw a nicely formatted exception on unrecoverable error
  [[noreturn]] void bail(const std::string& msg)
  {
    std::ostringstream oss;
    oss << "\n[FATAL] Fun4All_recoilJets :: " << msg << '\n';
    throw std::runtime_error(oss.str());
  }

  /// Trim whitespace from both ends (for robust list-file parsing)
  inline std::string trim(std::string s)
  {
    const char* ws = " \t\r\n";
    s.erase(0, s.find_first_not_of(ws));
    s.erase(s.find_last_not_of(ws) + 1);
    return s;
  }
}

namespace detail
{
  inline FastJetAlgoSub* fjAlgo(const float R)
  {
    FastJetOptions o{};               // IMPORTANT: value-initialize ALL fields to safe defaults
    o.algo            = Jet::ANTIKT;  // algorithm
    o.jet_R           = R;            // jet radius
    o.use_jet_min_pt  = true;         // enable a ptmin
    o.jet_min_pt      = 0.0f;         // ptmin value
    o.verbosity       = 0;            // quiet
    return new FastJetAlgoSub(o);
  }
}
//======================================================================
//  The actual steering macro
//======================================================================
void Fun4All_recoilJets(const int   nEvents   =  0,
                        const char* listFile  = "input_files.list",
                        const char* outRoot   = "TrigPlot.root",
                        const bool  verbose   = false)
{
  //--------------------------------------------------------------------
  // 0.  Banner & basic environment sanity
  //--------------------------------------------------------------------
  if (verbose) {
      std::cout << "\n>>> Fun4All_recoilJets – ana.495 driver <<<\n"
                << "    Input list : " << listFile  << '\n'
                << "    Output file: " << outRoot   << '\n'
                << "    nEvents    : " << nEvents   << (nEvents==0? " (all)\n":"\n");
  }

  Fun4AllServer* se = Fun4AllServer::instance();
  if (!se) detail::bail("unable to obtain Fun4AllServer instance!");

  //--------------------------------------------------------------------
  // 1.  Parse the file list & determine run / segment
  //--------------------------------------------------------------------
  std::ifstream list(listFile);
  if (!list.is_open())
      detail::bail("cannot open input list \"" + std::string(listFile) + "\"");

  std::vector<std::string> files;
  for (std::string line; std::getline(list, line); )
  {
      line = detail::trim(line);
      if (!line.empty()) files.emplace_back(line);
  }
  if (files.empty())
      detail::bail("input list \"" + std::string(listFile) + "\" is empty");

    const std::string& firstFile = files.front();

    // Simulation detection:
    //  - primary: RJ_DATASET=isSim
    //  - fallback: RJ_IS_SIM=1 (wrapper sets this)
    bool isSim = false;
    if (const char* ds = std::getenv("RJ_DATASET"))
    {
      std::string s = detail::trim(std::string(ds));
      std::transform(s.begin(), s.end(), s.begin(),
                     [](unsigned char c){ return std::tolower(c); });
      if (s == "issim" || s == "sim") isSim = true;
    }
    if (!isSim)
    {
      if (const char* f = std::getenv("RJ_IS_SIM"))
        isSim = (std::atoi(f) != 0);
    }

    auto runSegPair = Fun4AllUtils::GetRunSegment(firstFile);
    int run = runSegPair.first;
    int seg = runSegPair.second;

    if (run <= 0)
    {
      if (isSim)
      {
        // If the filename doesn't encode a run number, use a safe CDB timestamp.
        run = 1;
        seg = 0;

        // Optional override if you ever want it:
        //   export RJ_CDB_TIMESTAMP=XXXX
        if (const char* ts = std::getenv("RJ_CDB_TIMESTAMP"))
        {
          try
          {
            const int t = std::stoi(ts);
            if (t > 0) run = t;
          }
          catch (...) { /* keep default */ }
        }

        if (verbose)
          std::cout << "[INFO] Simulation mode: run/seg not in filename → using TIMESTAMP=" << run << "\n";
      }
      else
      {
        detail::bail("failed to extract run number from first file: " + firstFile);
      }
    }

    if (verbose)
          std::cout << "[INFO] Run=" << run << "  Seg=" << seg
                    << "  (" << files.size() << " files)\n";

  // --------------------------------------------------------------------
  // Global verbosity control (RJ_VERBOSITY from env; defaults to 10;
  // Condor detection → 0). Also silences std::cout/cerr globally when 0.
  // --------------------------------------------------------------------
  int vlevel = 10;
  if (const char* venv = std::getenv("RJ_VERBOSITY"))
  {
      try { vlevel = std::stoi(venv); } catch (...) { vlevel = 10; }
  }
  else
  {
      // Detect Condor jobs defensively; default to quiet in batch.
      if (std::getenv("_CONDOR_SCRATCH_DIR") || std::getenv("_CONDOR_JOB_AD"))
        vlevel = 0;
  }

  // RAII silence for global std::cout/cerr if vlevel==0
  struct ScopedSilence {
    std::ofstream   sink;
    std::streambuf* cout_save = nullptr;
    std::streambuf* cerr_save = nullptr;
    bool active = false;
    void enable() {
      if (active) return;
      sink.open("/dev/null");
      cout_save = std::cout.rdbuf(sink.rdbuf());
      cerr_save = std::cerr.rdbuf(sink.rdbuf());
      active = true;
    }
    ~ScopedSilence() {
      if (active) {
        std::cout.rdbuf(cout_save);
        std::cerr.rdbuf(cerr_save);
      }
    }
  } _silence;

  if (vlevel == 0) _silence.enable();

  // --------------------------------------------------------------------
  // 2.  CDB + IO managers
  // --------------------------------------------------------------------
  recoConsts* rc = recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG","newcdbtag");

  // Use run number as timestamp for real data.
  // For simulation, "run" like 28 is not a valid CDB timestamp for calo calibrations.
  unsigned long long cdbts = static_cast<unsigned long long>(run);

  if (isSim)
  {
      // Pick a known-good pp timestamp that has EMCal calibration payloads available.
      cdbts = 47289ULL;

      // Optional override:
      //   export RJ_CDB_TIMESTAMP=47289
      if (const char* ts = std::getenv("RJ_CDB_TIMESTAMP"))
      {
        char* end = nullptr;
        unsigned long long tmp = std::strtoull(ts, &end, 10);
        if (end != ts && tmp > 0ULL) cdbts = tmp;
      }
  }

  rc->set_uint64Flag("TIMESTAMP", cdbts);
  CDBInterface::instance()->Verbosity(vlevel);

  auto* flag = new FlagHandler();
  se->registerSubsystem(flag);

  auto* inDST = new Fun4AllDstInputManager("DSTcalofitting");
  for (const auto& f : files) inDST->AddFile(f);
  se->registerInputManager(inDST);

    // --------------------------------------------------------------------
    // 3.  Geometry + status + calibration + clustering
    // --------------------------------------------------------------------
    //
    // IMPORTANT:
    // RawClusterBuilderTemplate::process_event() ALWAYS requires "TOWERGEOM_CEMC".
    // But CaloGeomMapping with UseDetailedGeometry(true) publishes ONLY
    // "TOWERGEOM_CEMC_DETAILED" for CEMC.
    // So we ALWAYS create the legacy node (TOWERGEOM_CEMC), and optionally
    // also create the detailed node (TOWERGEOM_CEMC_DETAILED).
    //
    bool useDetailedCemcGeom = true;  // keep your current behavior by default
    if (const char* env = std::getenv("RJ_DETAILED_CEMC_GEOM"))
    {
      useDetailedCemcGeom = (std::atoi(env) != 0);
    }

    // Always publish the legacy/simple CEMC geometry node: TOWERGEOM_CEMC
    {
      auto* geomCemcLegacy = new CaloGeomMapping("Geom_CEMC");
      geomCemcLegacy->set_detector_name("CEMC");
      geomCemcLegacy->set_UseDetailedGeometry(false);
      se->registerSubsystem(geomCemcLegacy);
    }

    // Optionally publish the detailed CEMC geometry node: TOWERGEOM_CEMC_DETAILED
    if (useDetailedCemcGeom)
    {
      auto* geomCemcDetailed = new CaloGeomMapping("Geom_CEMC_DETAILED");
      geomCemcDetailed->set_detector_name("CEMC");
      geomCemcDetailed->set_UseDetailedGeometry(true);
      se->registerSubsystem(geomCemcDetailed);
    }

    // HCAL nodes (detailed not supported; CaloGeomMapping will fall back internally)
    for (const std::string& det : {"HCALIN","HCALOUT"})
    {
      auto* geom = new CaloGeomMapping(("Geom_" + det).c_str());
      geom->set_detector_name(det);
      geom->set_UseDetailedGeometry(true);
      se->registerSubsystem(geom);
    }


  if (vlevel > 0) std::cout << "status setters" << std::endl;
  CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
  statusEMC->set_detector_type(CaloTowerDefs::CEMC);
  statusEMC->set_time_cut(1);
  se->registerSubsystem(statusEMC);

  CaloTowerStatus *statusHCalIn = new CaloTowerStatus("HCALINSTATUS");
  statusHCalIn->set_detector_type(CaloTowerDefs::HCALIN);
  statusHCalIn->set_time_cut(2);
  se->registerSubsystem(statusHCalIn);

  CaloTowerStatus *statusHCALOUT = new CaloTowerStatus("HCALOUTSTATUS");
  statusHCALOUT->set_detector_type(CaloTowerDefs::HCALOUT);
  statusHCALOUT->set_time_cut(2);
  se->registerSubsystem(statusHCALOUT);

  if (vlevel > 0) std::cout << "Calibrating EMCal" << std::endl;
  CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
  calibEMC->set_detector_type(CaloTowerDefs::CEMC);
  se->registerSubsystem(calibEMC);

  if (vlevel > 0) std::cout << "Calibrating OHcal" << std::endl;
  CaloTowerCalib *calibOHCal = new CaloTowerCalib("HCALOUT");
  calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
  se->registerSubsystem(calibOHCal);

  if (vlevel > 0) std::cout << "Calibrating IHcal" << std::endl;
  CaloTowerCalib *calibIHCal = new CaloTowerCalib("HCALIN");
  calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
  se->registerSubsystem(calibIHCal);

  if (vlevel > 0) std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->set_threshold_energy(0.070);
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->set_UseTowerInfo(1);
  ClusterBuilder->set_UseAltZVertex(1);
  se->registerSubsystem(ClusterBuilder);


  if (vlevel > 0) std::cout << "Calibrating MBD" << std::endl;
  std::unique_ptr<MbdReco> mbdreco = std::make_unique<MbdReco>();
  se->registerSubsystem(mbdreco.release());

  if (!isSim)
    {
      if (vlevel > 0) std::cout << "Calibrating ZDC" << std::endl;
      auto* zdcreco = new ZdcReco();
      zdcreco->set_zdc1_cut(0.0);
      zdcreco->set_zdc2_cut(0.0);
      se->registerSubsystem(zdcreco);
    }
    else
    {
      if (vlevel > 0) std::cout << "[isSim] skipping ZdcReco (sim DST has no TOWERS_ZDC)" << std::endl;
  }


  if (vlevel > 0) std::cout << "Retrieving Vtx Info" << std::endl;
  std::unique_ptr<GlobalVertexReco> gvertex = std::make_unique<GlobalVertexReco>();
  se->registerSubsystem(gvertex.release());

    // Decide dataset early so we can gate Au+Au-only modules (centrality)
    // IMPORTANT: isSim should behave like pp-style reconstruction here.
    bool isAuAuData = true;

    if (const char* env = std::getenv("RJ_DATASET"))
    {
        std::string s = detail::trim(std::string(env));
        std::string sLower = s;
        std::transform(sLower.begin(), sLower.end(), sLower.begin(),
                       [](unsigned char c){ return std::tolower(c); });

        if (sLower == "issim" || sLower == "sim")
        {
          isSim = true;        // ensure consistency even if RJ_IS_SIM wasn't set
          isAuAuData = false;  // sim uses pp-style chain
        }
        else if (sLower == "ispp" || sLower == "pp")
        {
          isAuAuData = false;
        }
        else if (sLower == "isauau" || sLower == "auau" || sLower == "aa")
        {
          isAuAuData = true;
        }
        else
        {
          // Unknown token: fallback heuristic based on known pp runs
          isAuAuData = (run > 53864);
        }
    }
    else
    {
        // Fallback heuristic based on known pp runs: ≤ 53864 → p+p, > 53864 → Au+Au
        isAuAuData = (run > 53864);
    }

  if (isAuAuData)
  {
      std::cout << "building minbias classifier" << std::endl;
      auto* mb = new MinimumBiasClassifier();
      mb->Verbosity(0);
      mb->setOverwriteScale(
            "/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CentralityScale/42/6b/426bc1b56ba544201b0213766bee9478_cdb_centrality_scale_54912.root");
      se->registerSubsystem(mb);
      
      if (vlevel > 0) std::cout << "building centrality classifier (Au+Au)" << std::endl;
      auto* cent = new CentralityReco();
      cent->Verbosity(0);
      cent->setOverwriteScale(
              "/cvmfs/sphenix.sdcc.bnl.gov/calibrations/sphnxpro/cdb/CentralityScale/42/6b/426bc1b56ba544201b0213766bee9478_cdb_centrality_scale_54912.root");
      se->registerSubsystem(cent);
  }
  else
  {
        if (vlevel > 0) std::cout << "[pp dataset] skipping CentralityReco" << std::endl;
  }

  // --------------------------------------------------------------------
  // 4.  Jet reconstruction
  //   • Au+Au: HI-style retower + ρ subtraction (as before)
  //   • p+p  : simple raw Anti-kT jets on calibrated towers (no retower)
  // --------------------------------------------------------------------
  if (isAuAuData)
  {
      auto banner = [&](const std::string& m)
      {
        if (vlevel > 0) std::cout << "\n[BG-SUB] >>> " << m << std::endl;
      };
      const int vLvl = vlevel;

      // (i) Retower EMCal
      banner("(i)  Retowering EMCal to 0.025×0.025 towers");
      auto* rcemc = new RetowerCEMC("RetowerCEMC");
      rcemc->set_towerinfo(true);
      rcemc->set_frac_cut(0.50);
      rcemc->set_towerNodePrefix("TOWERINFO_CALIB");
      rcemc->Verbosity(0);
      se->registerSubsystem(rcemc);

      // Ensure DST/TOWER branch exists
      class EnsureTowerNode final : public SubsysReco
      {
       public:
        int InitRun(PHCompositeNode* top) override
        {
          PHNodeIterator it(top);
          auto* dst = dynamic_cast<PHCompositeNode*>(it.findFirst("PHCompositeNode","DST"));
          if (!dst) return Fun4AllReturnCodes::ABORTRUN;
          PHNodeIterator dstIt(dst);
          auto* tw = dynamic_cast<PHCompositeNode*>(dstIt.findFirst("PHCompositeNode","TOWER"));
          if (!tw)
          {
            tw = new PHCompositeNode("TOWER");
            dst->addNode(tw);
            if (std::getenv("RJ_VERBOSITY") && std::atoi(std::getenv("RJ_VERBOSITY")) > 0)
              std::cout << "[EnsureTowerNode]  created DST/TOWER branch\n";
          }
          return Fun4AllReturnCodes::EVENT_OK;
        }
      };
      se->registerSubsystem(new EnsureTowerNode());

      // (ii) RAW jets + (iii) backgrounds + (iv) subtract + (v) SUB1 jets + (vi) residual
      for (const auto& jnm : RecoilJets::kJetRadii)
      {
        const std::string radKey = jnm.key;           // dataset-aware map key (e.g., "r02")
        const int   D = std::stoi(radKey.substr(1));
        const float R = 0.1f * D;

        const std::string rawRecoName = "HIRecoSeedsRaw_" + radKey;
        const std::string subRecoName = "HIRecoSeedsSub_" + radKey;
        const std::string algoRaw     = "AntiKt_TowerInfo_" + rawRecoName;

        banner("(ii)  Reconstructing RAW Anti-kT jets  R=" + std::to_string(R));
        auto* seedReco = new JetReco(rawRecoName);
        seedReco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER, "TOWERINFO_CALIB"));
        seedReco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO,       "TOWERINFO_CALIB"));
        seedReco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO,      "TOWERINFO_CALIB"));
        seedReco->add_algo(detail::fjAlgo(R), algoRaw);
        seedReco->set_algo_node("ANTIKT");
        seedReco->set_input_node("TOWERINFO_CALIB");
        seedReco->Verbosity(0);
        se->registerSubsystem(seedReco);

        if (radKey == std::string("r02"))
        {
          banner("(iii-a)  UE density ρ from RAW jets (Sub1)");
          auto* dtb1 = new DetermineTowerBackground("DetTowerBkg_Sub1_r02");
          dtb1->SetBackgroundOutputName("TowerInfoBackground_Sub1");
          dtb1->SetSeedType(0);
          dtb1->SetSeedJetD(D);
          dtb1->set_towerNodePrefix("TOWERINFO_CALIB");
          dtb1->Verbosity(0);
          se->registerSubsystem(dtb1);

          banner("(iii-b)  UE density ρ for tower subtraction (Sub2)");
          auto* dtb2 = new DetermineTowerBackground("DetTowerBkg_Sub2_r02");
          dtb2->SetBackgroundOutputName("TowerInfoBackground_Sub2");
          dtb2->SetSeedType(0);
          dtb2->SetSeedJetD(D);
          dtb2->set_towerNodePrefix("TOWERINFO_CALIB");
          dtb2->Verbosity(0);
          se->registerSubsystem(dtb2);

          banner("(iv)  Subtracting ρ·A tower-by-tower  → *_SUB1");
          auto* st = new SubtractTowers("SubtractTowers_Sub1");
          st->set_towerinfo(true);
          st->set_towerNodePrefix("TOWERINFO_CALIB");
          st->Verbosity(0);
          se->registerSubsystem(st);
        }

        banner("(v)  Reconstructing jets on UE-subtracted towers  R=" + std::to_string(R));
        auto* subReco = new JetReco(subRecoName);
        subReco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_SUB1,   "TOWERINFO_CALIB"));
        subReco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO_SUB1, "TOWERINFO_CALIB"));
        subReco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO_SUB1,"TOWERINFO_CALIB"));
        subReco->add_algo(detail::fjAlgo(R), jnm.aa_node);   // dataset-aware Au+Au node name
        subReco->set_algo_node("ANTIKT");
        subReco->set_input_node("TOWERINFO_CALIB");
        subReco->Verbosity(0);
        se->registerSubsystem(subReco);

        if (radKey == std::string("r02"))
        {
          banner("(vi-a)  UE density ρ from SUB1 jets (Sub3)");
          auto* dtb3 = new DetermineTowerBackground("DetTowerBkg_Sub3_r02");
          dtb3->SetBackgroundOutputName("TowerInfoBackground_Sub3");
          dtb3->SetSeedType(1);
          dtb3->SetSeedJetD(D);
          dtb3->set_towerNodePrefix("TOWERINFO_CALIB");
          dtb3->Verbosity(0);
          se->registerSubsystem(dtb3);

          banner("(vi-b)  Copy jets + residual subtraction");
          auto* casj = new CopyAndSubtractJets("CopyAndSubtractJets_r02");
          casj->set_towerinfo(true);
          casj->set_towerNodePrefix("TOWERINFO_CALIB");
          casj->Verbosity(0);
          se->registerSubsystem(casj);
        }
      } // radii
  }   // Au+Au chain
  else
  {
    // p+p: simple raw Anti-kT jets on calibrated towers (no retower/ρ subtraction)
    auto banner = [&](const std::string& m)
    {
      if (vlevel > 0) std::cout << "\n[PP-JETS] >>> " << m << std::endl;
    };
    for (const auto& jnm : RecoilJets::kJetRadii)
    {
      const std::string radKey = jnm.key;
      const int   D = std::stoi(radKey.substr(1));
      const float R = 0.1f * D;

      const std::string recoName = "PPJetsRaw_" + radKey;

      banner("Reconstructing raw Anti-kT jets  R=" + std::to_string(R));
      auto* jpp = new JetReco(recoName);
      jpp->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO,    "TOWERINFO_CALIB"));
      jpp->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO,  "TOWERINFO_CALIB"));
      jpp->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO, "TOWERINFO_CALIB"));
      jpp->add_algo(detail::fjAlgo(R), jnm.pp_node);           // dataset-aware pp node name
      jpp->set_algo_node("ANTIKT");
      jpp->set_input_node("TOWERINFO_CALIB");
      jpp->Verbosity(0);
      se->registerSubsystem(jpp);
    }
  }

  // --------------------------------------------------------------------
  // 5.  Run-information helper (optional but handy)
  // --------------------------------------------------------------------
    if (!isSim)
    {
      auto* trigInfo = new TriggerRunInfoReco();
      trigInfo->Verbosity(vlevel);
      se->registerSubsystem(trigInfo);
    }
    else
    {
      if (vlevel > 0) std::cout << "[isSim] skipping TriggerRunInfoReco" << std::endl;
    }

  // Compute cluster isolation once per event (writes RawCluster::set_et_iso)
  // Only compute UE-subtracted isolation in Au+Au (SUB1 nodes exist there).
  const bool do_subtracted = isAuAuData;
  auto* clIso = new ClusterIso("ClusterIso",
                                 /*eTCut=*/0.0f,
                                 /*coneSize=*/3,           // R=0.3
                                 /*do_subtracted=*/do_subtracted,
                                 /*do_unsubtracted=*/true);
  clIso->Verbosity(0);
  se->registerSubsystem(clIso);
    
  // Build photon clusters
  auto* photonBuilder = new PhotonClusterBuilder("PhotonClusterBuilder");
  photonBuilder->set_input_cluster_node("CLUSTERINFO_CEMC");
  photonBuilder->set_output_photon_node("PHOTONCLUSTER_CEMC");
  photonBuilder->set_ET_threshold(0.07);
  photonBuilder->Verbosity(0);
  se->registerSubsystem(photonBuilder);

  auto* recoilJets = new RecoilJets(outRoot);
  recoilJets->setVzCut(30.);
  recoilJets->enableVzCut(true);
  recoilJets->setIsolationWP(1.08128, 0.0299107, 1.0,   0.30,  0.0);
  // Use the already-defined global vlevel
  recoilJets->setVerbose(vlevel);
  if (verbose) std::cout << "[INFO] RJ_VERBOSITY → " << vlevel << '\n';
  // Pick analysis type for the module (isPP / isAuAu / isSim), case-insensitive.
  // Fallback: ≤ 53864 → isPP, > 53864 → isAuAu.
  std::string dtype = "isAuAu";
  if (const char* env = std::getenv("RJ_DATASET"))
  {
            std::string s = detail::trim(std::string(env));
            std::string sLower = s;
            std::transform(sLower.begin(), sLower.end(), sLower.begin(),
                           [](unsigned char c){ return std::tolower(c); });

            if (sLower == "issim" || sLower == "sim")      dtype = "isSim";
            else if (sLower == "ispp" || sLower == "pp")   dtype = "isPP";
            else                                           dtype = "isAuAu";
    }
    else
    {
            dtype = (run <= 53864) ? "isPP" : "isAuAu";
  }

  if (verbose) std::cout << "[INFO] RJ_DATASET → " << dtype << '\n';
  recoilJets->setDataType(dtype);

  se->registerSubsystem(recoilJets);

  //--------------------------------------------------------------------
  // 6.  Run
  //--------------------------------------------------------------------
  try
  {
    if (verbose) std::cout << "[INFO] Starting event loop …\n";
    se->run(nEvents);
    se->End();
    if (verbose) std::cout << "[INFO] Finished successfully.\n";
  }
  catch (const std::exception& e)
  {
    detail::bail(std::string("exception in Fun4All: ") + e.what());
  }

//  //--------------------------------------------------------------------
//  // 7.  Clean exit
//  //--------------------------------------------------------------------
//  gSystem->Exit(0);
}

#endif   // ROOT_VERSION guard
