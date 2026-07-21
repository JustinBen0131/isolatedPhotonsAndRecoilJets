// Canary-only Arm B wrapper.  It preserves the authoritative unified steering
// and changes only the lifetime of neutral THE-106 observation registrations.
// The diagnostic writer implementation must be compiled as libTHE106C0R.so
// and made visible through LD_LIBRARY_PATH before this macro is loaded.

#pragma once

#define RJ_UNIFIED_ANALYSIS_AUAU 1
#include "../../../../macros/Fun4All_recoilJets_unified_impl.C"
#include "the106_c0r_diagnostic_writer.h"

#include <stdexcept>
#include <string>

R__LOAD_LIBRARY(libTHE106C0R.so)

void Fun4All_the106_c0r_AuAu(
    const int nEvents,
    const char* listFile,
    const char* outRoot,
    const char* cacheDirectory,
    const char* canaryId,
    const char* scientificContractSha256,
    const char* sourceManifestSha256,
    const char* configurationSha256,
    const char* reconstructionIdentity,
    const char* cdbIdentity,
    const char* modelReference,
    const char* modelSha256,
    const char* featureOrderSha256,
    const char* preprocessingIdentity,
    const char* runtimeProviderIdentity,
    const char* sourceRange,
    const bool verbose = false)
{
  the106::c0r::RunMetadata metadata;
  metadata.canary_id = canaryId ? canaryId : "";
  metadata.scientific_contract_sha256 =
      scientificContractSha256 ? scientificContractSha256 : "";
  metadata.source_manifest_sha256 =
      sourceManifestSha256 ? sourceManifestSha256 : "";
  metadata.configuration_sha256 =
      configurationSha256 ? configurationSha256 : "";
  metadata.reconstruction_identity =
      reconstructionIdentity ? reconstructionIdentity : "";
  metadata.cdb_identity = cdbIdentity ? cdbIdentity : "";
  metadata.model_reference = modelReference ? modelReference : "";
  metadata.model_sha256 = modelSha256 ? modelSha256 : "";
  metadata.feature_order_sha256 =
      featureOrderSha256 ? featureOrderSha256 : "";
  metadata.preprocessing_identity =
      preprocessingIdentity ? preprocessingIdentity : "";
  metadata.runtime_provider_identity =
      runtimeProviderIdentity ? runtimeProviderIdentity : "";
  metadata.source_range = sourceRange ? sourceRange : "";

  // Construction fails before Fun4All starts if the immutable cache namespace
  // is missing, nonempty, or not writable.  The registrations are thread-local
  // and default-null outside this lexical scope.
  the106::c0r::DiagnosticWriter writer(cacheDirectory ? cacheDirectory : "", metadata);
  if (!writer.good())
  {
    throw std::runtime_error("THE-106 C0-R diagnostic writer did not activate");
  }

  Fun4All_recoilJets_unified_impl(nEvents, listFile, outRoot, verbose);

  if (!writer.finish())
  {
    throw std::runtime_error(
        std::string("THE-106 C0-R diagnostic cache failed to seal: ") +
        writer.error());
  }
}
