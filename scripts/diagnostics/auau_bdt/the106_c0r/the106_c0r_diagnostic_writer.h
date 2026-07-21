#ifndef THE106_C0R_DIAGNOSTIC_WRITER_H
#define THE106_C0R_DIAGNOSTIC_WRITER_H

#include <phool/THE106Observation.h>

#include <cstdint>
#include <fstream>
#include <memory>
#include <string>

namespace the106
{
namespace c0r
{

// Canary-only provenance supplied by the sealed C0-R driver.  These values
// are recorded, never interpreted as scientific authority by the writer.
struct RunMetadata
{
  std::string canary_id;
  std::string scientific_contract_sha256;
  std::string source_manifest_sha256;
  std::string configuration_sha256;
  std::string reconstruction_identity;
  std::string cdb_identity;
  std::string model_reference;
  std::string model_sha256;
  std::string feature_order_sha256;
  std::string preprocessing_identity;
  std::string runtime_provider_identity;
  std::string source_range;
};

// A deliberately bounded observer/writer for the AU-AU-PHOTON-C0-RAW-QA
// canary.  It owns only diagnostic text files in an already-created empty
// directory and two thread-scoped observation registrations.  It never owns
// or mutates Fun4All, RecoilJets, ROOT output, input managers, or physics
// objects.
class DiagnosticWriter
{
 public:
  DiagnosticWriter(const std::string& cache_directory,
                   const RunMetadata& metadata);
  ~DiagnosticWriter() noexcept;

  DiagnosticWriter(const DiagnosticWriter&) = delete;
  DiagnosticWriter& operator=(const DiagnosticWriter&) = delete;

  bool good() const noexcept;
  const std::string& error() const noexcept;
  bool finish() noexcept;

  std::uint64_t serializationFailures() const noexcept;

 private:
  struct Streams;

  static void sourceCallback(void*, const c0h2::SourceObservation&) noexcept;
  static void scoreCallback(void*, const c0h2::ScoreObservation&) noexcept;
  static void eventCallback(void*, const c0rh::EventObservation&) noexcept;
  static void candidateCallback(void*, const c0rh::RawQACandidateObservation&) noexcept;
  static void fillCallback(void*, const c0rh::RawQAFillObservation&) noexcept;

  void writeSource(const c0h2::SourceObservation&);
  void writeScore(const c0h2::ScoreObservation&);
  void writeEvent(const c0rh::EventObservation&);
  void writeCandidate(const c0rh::RawQACandidateObservation&);
  void writeFill(const c0rh::RawQAFillObservation&);
  void fail(const char*) noexcept;

  std::string m_cache_directory;
  RunMetadata m_metadata;
  std::unique_ptr<Streams> m_streams;
  std::unique_ptr<c0h2::ScopedObservationRegistration> m_source_score_registration;
  std::unique_ptr<c0rh::ScopedObservationRegistration> m_recoil_registration;
  std::string m_error;
  std::uint64_t m_sequence;
  std::uint64_t m_score_sequence;
  std::uint64_t m_candidate_sequence;
  std::uint64_t m_serialization_failures;
  bool m_finished;
};

std::string percentEncode(const std::string& value);
std::string percentDecode(const std::string& value);
std::string binary32Hex(float value);
std::string binary64Hex(double value);
bool binary32FromHex(const std::string& text, float& value) noexcept;
bool binary64FromHex(const std::string& text, double& value) noexcept;

}  // namespace c0r
}  // namespace the106

#endif
