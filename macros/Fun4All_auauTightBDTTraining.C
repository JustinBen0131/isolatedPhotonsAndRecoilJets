#pragma once

#include <TSystem.h>

#define RJ_UNIFIED_ANALYSIS_AUAU 1
#include "Fun4All_recoilJets_unified_impl.C"

void Fun4All_auauTightBDTTraining(const int nEvents = 1000,
                                  const char* inputList = "",
                                  const char* outputFile = "auau_tight_bdt_training.root",
                                  const bool check = false)
{
  gSystem->Setenv("RJ_AUAU_BDT_EXTRACT_ONLY", "1");
  gSystem->Setenv("RJ_AUAU_BDT_TRAINING_TREE", "1");
  gSystem->Setenv("RJ_DISABLE_ID_FANOUT", "1");
  gSystem->Setenv("RJ_DISABLE_ISO_CONE_INTERNALIZATION", "1");
  gSystem->Setenv("RJ_DISABLE_JET_PT_INTERNALIZATION", "1");
  gSystem->Setenv("RJ_DISABLE_DPHI_INTERNALIZATION", "1");
  gSystem->Setenv("RJ_DIRECT_DST_DOALL", "1");
  gSystem->Setenv("RJ_INTERNAL_ISO_VIEWS", "");

  Fun4All_recoilJets_unified_impl(nEvents, inputList, outputFile, check);
}
