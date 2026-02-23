#ifndef CALO_CALIB_H
#define CALO_CALIB_H

#include <caloreco/CaloTowerCalib.h>
#include <caloreco/CaloTowerStatus.h>
#include <caloreco/RawClusterBuilderTemplate.h>
#include <caloreco/RawClusterDeadHotMask.h>
#include <caloreco/RawClusterPositionCorrection.h>

#include <ffamodules/CDBInterface.h>
#include <ffamodules/FlagHandler.h>

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>  // for Fun4AllServer

#include <phool/recoConsts.h>

#include <TSystem.h>  // for gSystem

R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libfun4allutils.so)

void Process_Calo_Calib()
{
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc = recoConsts::instance();

  // ------------------------------------------------------------
  // DEBUG: print CDB settings up front (this is where CEMCSTATUS pulls from)
  // ------------------------------------------------------------
  {
    std::string gtag = "(unset)";
    if (rc->FlagExist("CDB_GLOBALTAG"))
    {
      gtag = rc->get_StringFlag("CDB_GLOBALTAG");
    }
    std::cout << "[Calo_Calib][DBG] CDB_GLOBALTAG=" << gtag
              << "  TIMESTAMP=" << rc->get_uint64Flag("TIMESTAMP") << std::endl;

    int cdbV = 0;
    if (const char* env = std::getenv("RJ_CDB_VERBOSE")) cdbV = std::atoi(env);
    CDBInterface::instance()->Verbosity(cdbV);
    std::cout << "[Calo_Calib][DBG] CDBInterface Verbosity=" << cdbV << std::endl;
  }

  /////////////////
  // set MC or data
  bool isSim = true;
  int data_sim_runnumber_thres = 1000;
  if (rc->get_uint64Flag("TIMESTAMP") > data_sim_runnumber_thres)
  {
    isSim = false;
  }
  std::cout << "Calo Calib uses runnumber " << rc->get_uint64Flag("TIMESTAMP") << std::endl;

  //////////////////////
  // Input geometry node
  std::cout << "Adding Geometry file" << std::endl;
    Fun4AllInputManager *ingeo = new Fun4AllRunNodeInputManager("DST_GEO");
    std::string geoLocation = CDBInterface::instance()->getUrl("calo_geo");
    std::cout << "[Calo_Calib][DBG] calo_geo URL = '" << geoLocation << "'" << std::endl;
    if (geoLocation.empty())
    {
      std::cout << "[Calo_Calib][FATAL] CDBInterface::getUrl('calo_geo') returned empty string" << std::endl;
      gSystem->Exit(255);
    }
    ingeo->AddFile(geoLocation);
    se->registerInputManager(ingeo);

  //////////////////////////////
  // set statuses on raw towers
  std::cout << "status setters" << std::endl;
    // Optional: crank status-module verbosity (these are dying in InitRun)
    int statusV = 0;
    if (const char* env = std::getenv("RJ_STATUS_VERBOSE")) statusV = std::atoi(env);

    // DEBUG: node existence check before CEMCSTATUS InitRun
    class CaloStatusNodeProbe final : public SubsysReco
    {
     public:
      CaloStatusNodeProbe(const std::string& n) : SubsysReco(n) {}

      int InitRun(PHCompositeNode* topNode) override
      {
        printCheck("InitRun", topNode);
        return Fun4AllReturnCodes::EVENT_OK;
      }

      int process_event(PHCompositeNode* topNode) override
      {
        if (m_done) return Fun4AllReturnCodes::EVENT_OK;
        printCheck("process_event(evt1)", topNode);
        m_done = true;
        return Fun4AllReturnCodes::EVENT_OK;
      }

     private:
      void printCheck(const char* where, PHCompositeNode* topNode)
      {
        auto have = [&](const char* node) -> bool
        {
          return (findNode::getClass<TowerInfoContainer>(topNode, node) != nullptr);
        };

        std::cout << "[CaloStatusNodeProbe] " << where << " node check:"
                  << " TOWERINFO=" << (have("TOWERINFO") ? "OK" : "MISS")
                  << " TOWERINFO_CEMC=" << (have("TOWERINFO_CEMC") ? "OK" : "MISS")
                  << " TOWERINFO_CALIB=" << (have("TOWERINFO_CALIB") ? "OK" : "MISS")
                  << " TOWERINFO_CALIB_CEMC=" << (have("TOWERINFO_CALIB_CEMC") ? "OK" : "MISS")
                  << std::endl;
      }

      bool m_done = false;
    };

    se->registerSubsystem(new CaloStatusNodeProbe("CaloStatusNodeProbe_beforeCEMCSTATUS"));

    CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
    statusEMC->Verbosity(statusV);
    statusEMC->set_detector_type(CaloTowerDefs::CEMC);
    statusEMC->set_time_cut(1);

    std::cout << "[Calo_Calib][DBG] CEMCSTATUS configured:"
              << " isSim=" << (isSim ? "true" : "false")
              << " time_cut=1"
              << " Verbosity=" << statusV
              << std::endl;

    // MC Towers Status
    if(isSim) {
      // Uses threshold of 50% for towers be considered frequently bad.
      std::string calibName_hotMap = "CEMC_hotTowers_status";
      /* Systematic options (to be used as needed). */
      /* Uses threshold of 40% for towers be considered frequently bad. */
      // std::string calibName_hotMap = "CEMC_hotTowers_status_40";

      /* Uses threshold of 60% for towers be considered frequently bad. */
      // std::string calibName_hotMap = "CEMC_hotTowers_status_60";

      std::string calibdir = CDBInterface::instance()->getUrl(calibName_hotMap);
      std::cout << "[Calo_Calib][DBG] " << calibName_hotMap << " URL = '" << calibdir << "'" << std::endl;
      statusEMC->set_directURL_hotMap(calibdir);
    }
    se->registerSubsystem(statusEMC);

    CaloTowerStatus *statusHCalIn = new CaloTowerStatus("HCALINSTATUS");
    statusHCalIn->Verbosity(statusV);
    statusHCalIn->set_detector_type(CaloTowerDefs::HCALIN);
    statusHCalIn->set_time_cut(2);
    std::cout << "[Calo_Calib][DBG] HCALINSTATUS configured: time_cut=2 Verbosity=" << statusV << std::endl;
    se->registerSubsystem(statusHCalIn);

    CaloTowerStatus *statusHCALOUT = new CaloTowerStatus("HCALOUTSTATUS");
    statusHCALOUT->Verbosity(statusV);
    statusHCALOUT->set_detector_type(CaloTowerDefs::HCALOUT);
    statusHCALOUT->set_time_cut(2);
    std::cout << "[Calo_Calib][DBG] HCALOUTSTATUS configured: time_cut=2 Verbosity=" << statusV << std::endl;
    se->registerSubsystem(statusHCALOUT);

  ////////////////////
  // Calibrate towers
  std::cout << "Calibrating EMCal" << std::endl;
  CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
  calibEMC->set_detector_type(CaloTowerDefs::CEMC);
  se->registerSubsystem(calibEMC);

  std::cout << "Calibrating OHcal" << std::endl;
  CaloTowerCalib *calibOHCal = new CaloTowerCalib("HCALOUT");
  calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
  se->registerSubsystem(calibOHCal);

  std::cout << "Calibrating IHcal" << std::endl;
  CaloTowerCalib *calibIHCal = new CaloTowerCalib("HCALIN");
  calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
  se->registerSubsystem(calibIHCal);

  ////////////////
  // MC Calibration
  if (isSim && rc->get_uint64Flag("TIMESTAMP")<28) //in run28 and beyond we moved the MC calibration into the waveformsim module for data embedding
  {
    std::string MC_Calib = CDBInterface::instance()->getUrl("CEMC_MC_RECALIB");
    if (MC_Calib.empty())
    {
      std::cout << "No MC calibration found :( )" << std::endl;
      gSystem->Exit(0);
    }
    CaloTowerCalib *calibEMC_MC = new CaloTowerCalib("CEMCCALIB_MC");
    calibEMC_MC->set_detector_type(CaloTowerDefs::CEMC);
    calibEMC_MC->set_inputNodePrefix("TOWERINFO_CALIB_");
    calibEMC_MC->set_outputNodePrefix("TOWERINFO_CALIB_");
    calibEMC_MC->set_directURL(MC_Calib);
    calibEMC_MC->set_doCalibOnly(true);
    se->registerSubsystem(calibEMC_MC);
  }

  //////////////////
  // Clusters
  std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->set_threshold_energy(0.070);  // for when using basic calibration
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
  ClusterBuilder->set_UseAltZVertex(1); // Use MBD Vertex for vertex-based corrections
  se->registerSubsystem(ClusterBuilder);

  // currently NOT included!
  // std::cout << "Applying Position Dependent Correction" << std::endl;
  // RawClusterPositionCorrection *clusterCorrection = new RawClusterPositionCorrection("CEMC");
  // clusterCorrection->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
  // se->registerSubsystem(clusterCorrection);
}

#endif
