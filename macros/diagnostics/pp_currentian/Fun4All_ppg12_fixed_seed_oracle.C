#ifndef FUN4ALL_PPG12_FIXED_SEED_ORACLE_C
#define FUN4ALL_PPG12_FIXED_SEED_ORACLE_C

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <TROOT.h>
#include <TSystem.h>

#include <cstdlib>
#include <array>
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
constexpr std::array<unsigned int, 5> kHistoricalSeeds = {
    2991264730U,
    4256268992U,
    2394322166U,
    874466025U,
    2240380304U};
constexpr unsigned int kHistoricalPedestalSequence = 534U;

void load_historical_seed_contract()
{
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED"))
  {
    throw std::runtime_error(
        "historical PPG12 replay forbids recoConsts RANDOMSEED");
  }
  PHRandomSeed::Verbosity(1);
  for (const unsigned int seed : kHistoricalSeeds)
  {
    PHRandomSeed::LoadSeed(seed);
  }
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
    const std::string &g4List = "g4hits_first5.list",
    const std::string &truthJetList = "dst_truth_jet_first5.list",
    const std::string &outputFile = "output_sim.root")
{
  ppg12_paired_oracle::load_historical_seed_contract();
  ppg12_paired_oracle::emit_runtime_fingerprint("ppg12");
  std::cout
      << "ORACLE_SEED_CONTRACT side=ppg12 mode=historical_fifo_replay_v2"
      << " rc_randomseed=absent"
      << " ph_seed_sequence=2991264730,4256268992,2394322166,874466025,2240380304"
      << " pedestal_sequence="
      << ppg12_paired_oracle::kHistoricalPedestalSequence << std::endl;
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
