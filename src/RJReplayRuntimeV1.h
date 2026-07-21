#ifndef RJ_REPLAY_RUNTIME_V1_H
#define RJ_REPLAY_RUNTIME_V1_H

#include "RJReplayFoundationV1.h"

#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

namespace RJReplayRuntimeV1
{
using namespace RJReplayFoundationV1;

inline std::string env(const char* key)
{
  const char* value = std::getenv(key);
  return value ? std::string(value) : std::string{};
}

inline bool envEnabled(const char* key)
{
  const std::string value = env(key);
  return value == "1" || value == "true" || value == "TRUE" || value == "yes" || value == "on";
}

inline int envInt(const char* key, int fallback)
{
  const std::string value = env(key);
  if (value.empty()) return fallback;
  try { return std::stoi(value); }
  catch (...) { return fallback; }
}

struct EventBundle
{
  EventRow event;
  std::vector<PhotonCandidateRow> candidates;
  std::vector<ModelEvaluationRow> models;
  std::vector<ShowerCellRow> shower_cells;
  std::vector<IsolationConstituentRow> isolation_constituents;
  std::vector<IsolationWitnessRow> isolation_witnesses;
  std::vector<JetRow> jets;
  std::vector<JetConstituentRow> jet_constituents;
  std::vector<PhotonJetPairRow> pairs;
  std::vector<TruthPhotonRow> truth_photons;
  std::vector<TruthJetRow> truth_jets;
  std::vector<RecoTruthLinkRow> links;
  std::vector<WeightComponentRow> weights;
  std::vector<EventDisplaySnapshotRow> snapshots;
};

class Runtime
{
 public:
  bool initialize(TFile* file, const SourceOccurrenceRow& source, std::string* error = nullptr)
  {
    Metadata metadata;
    metadata.schema_sha256 = env("RJ_REPLAY_SCHEMA_SHA256");
    metadata.semantic_sha256 = env("RJ_REPLAY_SEMANTIC_SHA256");
    metadata.source_sha256 = env("RJ_REPLAY_SOURCE_SHA256");
    metadata.model_sha256 = env("RJ_REPLAY_MODEL_SHA256");
    metadata.config_sha256 = env("RJ_REPLAY_CONFIG_SHA256");
    metadata.code_sha256 = env("RJ_REPLAY_CODE_SHA256");
    if (!m_writer.initialize(file, metadata, error)) return false;
    if (!m_writer.fill(source, error)) return false;
    m_source_id = source.id;
    return true;
  }

  bool write(EventBundle& bundle, std::string* error = nullptr)
  {
    if (bundle.event.source_id.isNull()) bundle.event.source_id = m_source_id;
    if (!m_writer.fill(bundle.event, error)) return false;
    for (const auto& row : bundle.candidates) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.models) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.shower_cells) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.isolation_constituents) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.isolation_witnesses) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.jets) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.jet_constituents) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.pairs) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.truth_photons) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.truth_jets) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.links) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.weights) if (!m_writer.fill(row, error)) return false;
    for (const auto& row : bundle.snapshots) if (!m_writer.fill(row, error)) return false;
    return true;
  }

  bool finish(std::string* error = nullptr) { return m_writer.finish(error); }

 private:
  Writer m_writer;
  Identity128 m_source_id;
};

template <class Function>
class ScopeExit
{
 public:
  explicit ScopeExit(Function&& function) : m_function(std::move(function)) {}
  ScopeExit(const ScopeExit&) = delete;
  ScopeExit& operator=(const ScopeExit&) = delete;
  ~ScopeExit() { m_function(); }

 private:
  Function m_function;
};

template <class Function>
ScopeExit<Function> onScopeExit(Function&& function)
{
  return ScopeExit<Function>(std::forward<Function>(function));
}
}  // namespace RJReplayRuntimeV1

#endif
