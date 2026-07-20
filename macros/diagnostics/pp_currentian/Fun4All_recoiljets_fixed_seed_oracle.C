#ifndef FUN4ALL_RECOILJETS_FIXED_SEED_ORACLE_C
#define FUN4ALL_RECOILJETS_FIXED_SEED_ORACLE_C

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <TROOT.h>
#include <TSystem.h>

#include <cstdlib>
#include <iostream>
#include <stdexcept>

// The isolated, new.17-built RecoilJets macro is loaded by the worker first.
void Fun4All_recoilJets(
    const int nEvents,
    const char *listFile,
    const char *outRoot,
    const bool verbose);

namespace recoiljets_paired_oracle
{
void require_seed_contract(const int seed)
{
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED") && rc->get_IntFlag("RANDOMSEED") != seed)
  {
    throw std::runtime_error("conflicting recoConsts RANDOMSEED in RecoilJets oracle process");
  }
  rc->set_IntFlag("RANDOMSEED", seed);
  PHRandomSeed::Verbosity(1);
}

void emit_runtime_fingerprint(const char *side)
{
  const char *offline = std::getenv("OFFLINE_MAIN");
  const char *profile = std::getenv("RJ_PPG12_ORACLE_RUNTIME_PROFILE");
  const char *caloCalib = gSystem->Which(gROOT->GetMacroPath(), "Calo_Calib.C");

  std::cout << "ORACLE_RUNTIME side=" << side
            << " profile=" << (profile ? profile : "<unset>")
            << " offline_main=" << (offline ? offline : "<unset>") << std::endl;
  std::cout << "ORACLE_CALO_CALIB side=" << side
            << " path=" << (caloCalib ? caloCalib : "<missing>") << std::endl;
  for (const char *name : {"libcalo_reco.so", "libRecoilJets.so", "libclusteriso.so", "libjetbase.so"})
  {
    const char *path = gSystem->DynamicPathName(name, true);
    std::cout << "ORACLE_DYNAMIC_LIBRARY side=" << side << " name=" << name
              << " path=" << (path ? path : "<missing>") << std::endl;
  }
  std::cout << "ORACLE_LIBRARIES_BEGIN side=" << side << std::endl;
  std::cout << gSystem->GetLibraries() << std::endl;
  std::cout << "ORACLE_LIBRARIES_END side=" << side << std::endl;
}
}  // namespace recoiljets_paired_oracle

void Fun4All_recoiljets_fixed_seed_oracle(
    const int seed = 42,
    const char *combinedList = "recoil_first5.list",
    const char *outputFile = "recoil.root")
{
  recoiljets_paired_oracle::require_seed_contract(seed);
  recoiljets_paired_oracle::emit_runtime_fingerprint("recoiljets");
  std::cout << "ORACLE_SEED_CONTRACT side=recoiljets rc_seed=" << seed << std::endl;
  std::cout << "ORACLE_SOURCE_GRAPH side=recoiljets columns=NONE,g4,truthjet,NONE,NONE"
            << std::endl;
  Fun4All_recoilJets(0, combinedList, outputFile, true);
}

#endif
