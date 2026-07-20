#ifndef FUN4ALL_PPG12_FIXED_SEED_ORACLE_C
#define FUN4ALL_PPG12_FIXED_SEED_ORACLE_C

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <TROOT.h>
#include <TSystem.h>

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

// The preserved macro is loaded by the worker before this wrapper.
void Fun4All_run_sim(
    const int nEvents,
    const std::string &inputFile0,
    const std::string &inputFile1,
    const std::string &inputFile3,
    const std::string &inputFile4,
    const std::string &inputFile2,
    const std::string &outputFile,
    const std::string &outputDSTFile,
    const std::string &outDSTdir,
    const std::string &cdbtag);

namespace ppg12_paired_oracle
{
void require_seed_contract(const int seed)
{
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED") && rc->get_IntFlag("RANDOMSEED") != seed)
  {
    throw std::runtime_error("conflicting recoConsts RANDOMSEED in PPG12 oracle process");
  }
  rc->set_IntFlag("RANDOMSEED", seed);
  PHRandomSeed::Verbosity(1);
}

void emit_runtime_fingerprint(const char *side)
{
  const char *offline = std::getenv("OFFLINE_MAIN");
  const char *profile = std::getenv("RJ_PPG12_ORACLE_RUNTIME_PROFILE");
  const char *caloCalib = gSystem->Which(gROOT->GetMacroPath(), "Calo_Calib.C");
  const char *caloAna = gSystem->DynamicPathName("libCaloAna24.so", true);

  std::cout << "ORACLE_RUNTIME side=" << side
            << " profile=" << (profile ? profile : "<unset>")
            << " offline_main=" << (offline ? offline : "<unset>") << std::endl;
  std::cout << "ORACLE_CALO_CALIB side=" << side
            << " path=" << (caloCalib ? caloCalib : "<missing>") << std::endl;
  std::cout << "ORACLE_DYNAMIC_LIBRARY side=" << side
            << " name=libCaloAna24.so path=" << (caloAna ? caloAna : "<missing>")
            << std::endl;
  std::cout << "ORACLE_LIBRARIES_BEGIN side=" << side << std::endl;
  std::cout << gSystem->GetLibraries() << std::endl;
  std::cout << "ORACLE_LIBRARIES_END side=" << side << std::endl;
}
}  // namespace ppg12_paired_oracle

void Fun4All_ppg12_fixed_seed_oracle(
    const int seed = 42,
    const std::string &g4List = "g4hits_first5.list",
    const std::string &truthJetList = "dst_truth_jet_first5.list",
    const std::string &outputFile = "output_sim.root")
{
  ppg12_paired_oracle::require_seed_contract(seed);
  ppg12_paired_oracle::emit_runtime_fingerprint("ppg12");
  std::cout << "ORACLE_SEED_CONTRACT side=ppg12 rc_seed=" << seed << std::endl;
  std::cout << "ORACLE_SOURCE_GRAPH side=ppg12 columns=NONE,g4,truthjet,NONE,NONE"
            << std::endl;

  // Direct remote inspection of the frozen Photon5 SI macro establishes that
  // only INPUTREADHITS indices 0 and 4 are active.  The other arguments are
  // deliberately inert sentinels; passing live auxiliary lists here would
  // misrepresent the executable oracle.
  Fun4All_run_sim(
      0,
      g4List,
      "NONE",
      "NONE",
      truthJetList,
      "NONE",
      outputFile,
      "NONE.root",
      ".",
      "MDC2");
}

#endif
