#include "the106_c0r_diagnostic_writer.h"

#include <cerrno>
#include <cctype>
#include <cstring>
#include <fcntl.h>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unistd.h>

namespace the106
{
namespace c0r
{
namespace
{

const char* nullable(const char* value) noexcept
{
  return value ? value : "";
}

const char* boolean(bool value) noexcept
{
  return value ? "1" : "0";
}

std::string pointerHex(const void* value)
{
  std::ostringstream out;
  out << "0x" << std::hex << std::setfill('0')
      << std::setw(static_cast<int>(2 * sizeof(std::uintptr_t)))
      << reinterpret_cast<std::uintptr_t>(value);
  return out.str();
}

std::string unsignedHex(std::uint64_t value, std::size_t digits)
{
  std::ostringstream out;
  out << "0x" << std::hex << std::setfill('0')
      << std::setw(static_cast<int>(digits)) << value;
  return out.str();
}

bool parseHex(const std::string& text, std::uint64_t& value,
              std::size_t max_digits) noexcept
{
  try
  {
    if (text.size() < 3 || text.size() > max_digits + 2 ||
        text[0] != '0' || (text[1] != 'x' && text[1] != 'X'))
    {
      return false;
    }
    std::size_t used = 0;
    value = std::stoull(text.substr(2), &used, 16);
    return used == text.size() - 2;
  }
  catch (...)
  {
    return false;
  }
}

bool isSha256(const std::string& value) noexcept
{
  if (value.size() != 64) return false;
  for (std::string::const_iterator iter = value.begin(); iter != value.end(); ++iter)
    if (!std::isxdigit(static_cast<unsigned char>(*iter))) return false;
  return true;
}

void reserveExclusive(const std::string& path)
{
  const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0444);
  if (fd < 0)
  {
    std::ostringstream message;
    message << "refusing to overwrite diagnostic cache file " << path
            << ": " << std::strerror(errno);
    throw std::runtime_error(message.str());
  }
  if (::close(fd) != 0)
  {
    throw std::runtime_error("unable to close newly reserved cache file " + path);
  }
  // The files are private to an immutable canary namespace; make them
  // writable only for this process after the O_EXCL reservation.
  if (::chmod(path.c_str(), 0644) != 0)
  {
    throw std::runtime_error("unable to chmod newly reserved cache file " + path);
  }
}

std::ofstream openExclusive(const std::string& directory,
                            const std::string& name)
{
  const std::string path = directory + "/" + name;
  reserveExclusive(path);
  std::ofstream stream(path.c_str(), std::ios::out | std::ios::trunc);
  if (!stream)
  {
    throw std::runtime_error("unable to open diagnostic cache file " + path);
  }
  return stream;
}

void writeContext(std::ostream& out, const c0h2::CandidateContext& context)
{
  out << context.pair.generation << '\t'
      << context.pair.pair_ordinal << '\t'
      << context.delivered_event_ordinal << '\t'
      << boolean(context.run_valid) << '\t'
      << context.run_number << '\t'
      << boolean(context.event_valid) << '\t'
      << context.event_number << '\t'
      << context.container_key << '\t'
      << context.cluster_id << '\t'
      << context.producer_encounter_ordinal << '\t'
      << boolean(context.canonical_view) << '\t'
      << percentEncode(nullable(context.view_key)) << '\t'
      << percentEncode(nullable(context.module_name));
}

void writeOccurrence(std::ostream& out,
                     const c0h2::OccurrenceToken& occurrence)
{
  out << occurrence.generation << '\t'
      << occurrence.manager_ordinal << '\t'
      << occurrence.occurrence_ordinal;
}

}  // namespace

struct DiagnosticWriter::Streams
{
  explicit Streams(const std::string& directory)
    : metadata(openExclusive(directory, "metadata.tsv"))
    , source(openExclusive(directory, "source_observations.tsv"))
    , scores(openExclusive(directory, "score_observations.tsv"))
    , score_features(openExclusive(directory, "score_features.tsv"))
    , events(openExclusive(directory, "event_observations.tsv"))
    , event_triggers(openExclusive(directory, "event_active_triggers.tsv"))
    , candidates(openExclusive(directory, "rawqa_candidates.tsv"))
    , candidate_values(openExclusive(directory, "rawqa_candidate_values.tsv"))
    , candidate_triggers(openExclusive(directory, "rawqa_candidate_triggers.tsv"))
    , fills(openExclusive(directory, "rawqa_fill_witnesses.tsv"))
    , completion(openExclusive(directory, "completion.tsv"))
  {
    metadata << "key\tvalue\n";
    source << "sequence\tkind\toccurrence_generation\tmanager_ordinal"
              "\toccurrence_ordinal\tpair_generation\tpair_ordinal"
              "\tmanager_owner\tsync_owner\tio_owner\tmanager_name"
              "\tsource_descriptor\tsource_entry\tread_attempt_ordinal"
              "\tread_kind\tread_outcome\tphase\tterminal\ttransition"
              "\tnative_status\treason\n";
    scores << "sequence\tpair_generation\tpair_ordinal\tdelivered_event_ordinal"
              "\trun_valid\trun_number\tevent_valid\tevent_number"
              "\tcontainer_key\tcluster_id\tproducer_encounter_ordinal"
              "\tcanonical_view\tview_key\tmodule_name\tmodel_reference"
              "\tscoring_mode\tmodel_sha256\tfeature_order_sha256"
              "\tpreprocessing_identity\truntime_provider_identity\tstatus"
              "\tscore_valid\tscore_bits\tfailing_feature_index\tfeature_count\n";
    score_features << "score_sequence\tfeature_index\tfeature_name\tvalue_bits\n";
    events << "sequence\tkind\tpair_generation\tpair_ordinal\toccurrence_count"
              "\toccurrence0_generation\toccurrence0_manager_ordinal"
              "\toccurrence0_ordinal\toccurrence1_generation"
              "\toccurrence1_manager_ordinal\toccurrence1_ordinal"
              "\tdelivered_event_ordinal\trun_valid\trun_number"
              "\tevent_valid\tevent_number\tcentrality_valid\tcentrality_bits"
              "\toutcome\tnative_status\toutcome_reason"
              "\tencountered_candidate_count\tadmitted_candidate_count"
              "\tscored_candidate_count\trawqa_contributor_count"
              "\trawqa_fill_count\tactive_trigger_count\tzero_candidate\tzero_fill\tfinalized"
              "\tmodule_name\n";
    event_triggers << "event_sequence\ttrigger_index\ttrigger_name\n";
    candidates << "sequence\tpair_generation\tpair_ordinal\tdelivered_event_ordinal"
                  "\trun_valid\trun_number\tevent_valid\tevent_number"
                  "\tcontainer_key\tcluster_id\tproducer_encounter_ordinal"
                  "\tcanonical_view_context\tview_key\tmodule_name"
                  "\tphoton_pt_slice\tcentrality_slice\tcanonical_view"
                  "\tview_suffix\ttight_tag\tpreselection_pass"
                  "\tdirect_candidate_admitted\tvalid_mask\tvalue_count"
                  "\tactive_trigger_count\n";
    candidate_values << "candidate_sequence\tvalue_index\tvariable_name"
                        "\tvalue_bits\tvalid\n";
    candidate_triggers << "candidate_sequence\ttrigger_index\ttrigger_name\n";
    fills << "sequence\tpair_generation\tpair_ordinal\tdelivered_event_ordinal"
             "\trun_valid\trun_number\tevent_valid\tevent_number"
             "\tcontainer_key\tcluster_id\tproducer_encounter_ordinal"
             "\tcanonical_view_context\tview_key\tmodule_name\tvariable_name"
             "\ttrigger_name\ttag_name\tview_suffix\tdirectory_name"
             "\tobject_name\tobject_class\tobject_contract_key"
             "\tphoton_pt_slice\tcentrality_slice\tcanonical_view\ttight_tag"
             "\tvalue_bits\tweight_bits\tfill_ordinal\tvalue_valid\tfilled"
             "\tsumw2_enabled\n";
    completion << "key\tvalue\n";
  }

  bool good() const
  {
    return metadata.good() && source.good() && scores.good() &&
           score_features.good() && events.good() && event_triggers.good() && candidates.good() &&
           candidate_values.good() && candidate_triggers.good() &&
           fills.good() && completion.good();
  }

  void flush()
  {
    metadata.flush();
    source.flush();
    scores.flush();
    score_features.flush();
    events.flush();
    event_triggers.flush();
    candidates.flush();
    candidate_values.flush();
    candidate_triggers.flush();
    fills.flush();
    completion.flush();
  }

  std::ofstream metadata;
  std::ofstream source;
  std::ofstream scores;
  std::ofstream score_features;
  std::ofstream events;
  std::ofstream event_triggers;
  std::ofstream candidates;
  std::ofstream candidate_values;
  std::ofstream candidate_triggers;
  std::ofstream fills;
  std::ofstream completion;
};

std::string percentEncode(const std::string& value)
{
  static const char digits[] = "0123456789ABCDEF";
  std::string encoded;
  encoded.reserve(value.size());
  for (std::string::const_iterator iter = value.begin(); iter != value.end(); ++iter)
  {
    const unsigned char ch = static_cast<unsigned char>(*iter);
    const bool safe = std::isalnum(ch) || ch == '.' || ch == '_' ||
                      ch == '-' || ch == '/';
    if (safe)
    {
      encoded.push_back(static_cast<char>(ch));
    }
    else
    {
      encoded.push_back('%');
      encoded.push_back(digits[(ch >> 4) & 0x0f]);
      encoded.push_back(digits[ch & 0x0f]);
    }
  }
  return encoded;
}

std::string percentDecode(const std::string& value)
{
  std::string decoded;
  decoded.reserve(value.size());
  for (std::size_t index = 0; index < value.size(); ++index)
  {
    if (value[index] != '%')
    {
      decoded.push_back(value[index]);
      continue;
    }
    if (index + 2 >= value.size())
    {
      throw std::runtime_error("truncated percent escape");
    }
    const std::string digits = value.substr(index + 1, 2);
    std::size_t used = 0;
    const unsigned long byte = std::stoul(digits, &used, 16);
    if (used != 2 || byte > 255)
    {
      throw std::runtime_error("invalid percent escape");
    }
    decoded.push_back(static_cast<char>(byte));
    index += 2;
  }
  return decoded;
}

std::string binary32Hex(float value)
{
  std::uint32_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value), "IEEE-754 binary32 required");
  std::memcpy(&bits, &value, sizeof(bits));
  return unsignedHex(bits, 8);
}

std::string binary64Hex(double value)
{
  std::uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value), "IEEE-754 binary64 required");
  std::memcpy(&bits, &value, sizeof(bits));
  return unsignedHex(bits, 16);
}

bool binary32FromHex(const std::string& text, float& value) noexcept
{
  std::uint64_t parsed = 0;
  if (!parseHex(text, parsed, 8) || parsed > 0xffffffffULL) return false;
  const std::uint32_t bits = static_cast<std::uint32_t>(parsed);
  std::memcpy(&value, &bits, sizeof(value));
  return true;
}

bool binary64FromHex(const std::string& text, double& value) noexcept
{
  std::uint64_t bits = 0;
  if (!parseHex(text, bits, 16)) return false;
  std::memcpy(&value, &bits, sizeof(value));
  return true;
}

DiagnosticWriter::DiagnosticWriter(const std::string& cache_directory,
                                   const RunMetadata& metadata)
  : m_cache_directory(cache_directory)
  , m_metadata(metadata)
  , m_streams()
  , m_source_score_registration()
  , m_recoil_registration()
  , m_error()
  , m_sequence(0)
  , m_score_sequence(0)
  , m_candidate_sequence(0)
  , m_serialization_failures(0)
  , m_finished(false)
{
  struct stat state;
  if (::stat(cache_directory.c_str(), &state) != 0 || !S_ISDIR(state.st_mode))
  {
    throw std::runtime_error(
        "diagnostic cache directory must already exist: " + cache_directory);
  }
  if (metadata.canary_id.empty() ||
      !isSha256(metadata.scientific_contract_sha256) ||
      !isSha256(metadata.source_manifest_sha256) ||
      !isSha256(metadata.configuration_sha256) ||
      metadata.reconstruction_identity.empty() || metadata.cdb_identity.empty() ||
      metadata.model_reference.empty() || !isSha256(metadata.model_sha256) ||
      !isSha256(metadata.feature_order_sha256) ||
      metadata.preprocessing_identity.empty() ||
      metadata.runtime_provider_identity.empty() || metadata.source_range.empty())
  {
    throw std::runtime_error("incomplete or malformed sealed C0-R metadata");
  }

  m_streams.reset(new Streams(cache_directory));
  const std::pair<const char*, const std::string*> metadata_rows[] = {
      std::make_pair("format", static_cast<const std::string*>(nullptr)),
      std::make_pair("canary_id", &m_metadata.canary_id),
      std::make_pair("scientific_contract_sha256", &m_metadata.scientific_contract_sha256),
      std::make_pair("source_manifest_sha256", &m_metadata.source_manifest_sha256),
      std::make_pair("configuration_sha256", &m_metadata.configuration_sha256),
      std::make_pair("reconstruction_identity", &m_metadata.reconstruction_identity),
      std::make_pair("cdb_identity", &m_metadata.cdb_identity),
      std::make_pair("model_reference", &m_metadata.model_reference),
      std::make_pair("model_sha256", &m_metadata.model_sha256),
      std::make_pair("feature_order_sha256", &m_metadata.feature_order_sha256),
      std::make_pair("preprocessing_identity", &m_metadata.preprocessing_identity),
      std::make_pair("runtime_provider_identity", &m_metadata.runtime_provider_identity),
      std::make_pair("source_range", &m_metadata.source_range)};
  for (std::size_t index = 0;
       index < sizeof(metadata_rows) / sizeof(metadata_rows[0]); ++index)
  {
    const std::string value = metadata_rows[index].second
        ? *metadata_rows[index].second : "THE106_C0R_CACHE_V1";
    m_streams->metadata << metadata_rows[index].first << '\t'
                        << percentEncode(value) << '\n';
  }

  c0h2::ScoreRegistrationMetadata score_metadata;
  score_metadata.model_sha256 = metadata.model_sha256;
  score_metadata.feature_order_sha256 = metadata.feature_order_sha256;
  score_metadata.preprocessing_identity = metadata.preprocessing_identity;
  score_metadata.runtime_provider_identity = metadata.runtime_provider_identity;

  // c0h2 is registered first and therefore destroyed last.  This preserves
  // pair/delivery context through every c0rh event finalization callback.
  m_source_score_registration.reset(new c0h2::ScopedObservationRegistration(
      &DiagnosticWriter::sourceCallback, this,
      &DiagnosticWriter::scoreCallback, this,
      score_metadata));
  m_recoil_registration.reset(new c0rh::ScopedObservationRegistration(
      &DiagnosticWriter::eventCallback, this,
      &DiagnosticWriter::candidateCallback, this,
      &DiagnosticWriter::fillCallback, this));

  if (!m_source_score_registration->active() ||
      !m_recoil_registration->active() || !m_streams->good())
  {
    throw std::runtime_error("unable to activate complete C0-R observation writer");
  }
}

DiagnosticWriter::~DiagnosticWriter() noexcept
{
  finish();
}

bool DiagnosticWriter::good() const noexcept
{
  return m_error.empty() && m_streams && m_streams->good() &&
         !m_finished;
}

const std::string& DiagnosticWriter::error() const noexcept
{
  return m_error;
}

std::uint64_t DiagnosticWriter::serializationFailures() const noexcept
{
  return m_serialization_failures;
}

void DiagnosticWriter::fail(const char* message) noexcept
{
  ++m_serialization_failures;
  if (m_error.empty())
  {
    try
    {
      m_error = message ? message : "diagnostic serialization failure";
    }
    catch (...)
    {
    }
  }
}

bool DiagnosticWriter::finish() noexcept
{
  if (m_finished) return m_error.empty();
  m_finished = true;
  try
  {
    const c0h2::RegistrationStats source_stats =
        m_source_score_registration
            ? m_source_score_registration->stats() : c0h2::RegistrationStats();
    const c0rh::RegistrationStats recoil_stats =
        m_recoil_registration
            ? m_recoil_registration->stats() : c0rh::RegistrationStats();
    const bool closed_world = c0h2::validateObservationStateClosedWorld();

    m_recoil_registration.reset();
    m_source_score_registration.reset();

    m_streams->completion << "source_observer_failures\t"
                          << source_stats.observer_failures << '\n'
                          << "source_invariant_violations\t"
                          << source_stats.invariant_violations << '\n'
                          << "source_callbacks\t" << source_stats.source_callbacks << '\n'
                          << "score_callbacks\t" << source_stats.score_callbacks << '\n'
                          << "recoil_observer_failures\t"
                          << recoil_stats.observer_failures << '\n'
                          << "recoil_invariant_violations\t"
                          << recoil_stats.invariant_violations << '\n'
                          << "event_callbacks\t" << recoil_stats.event_callbacks << '\n'
                          << "candidate_callbacks\t" << recoil_stats.candidate_callbacks << '\n'
                          << "fill_callbacks\t" << recoil_stats.fill_callbacks << '\n'
                          << "closed_world\t" << boolean(closed_world) << '\n'
                          << "serialization_failures\t"
                          << m_serialization_failures << '\n';
    m_streams->flush();
    if (!closed_world)
    {
      fail("source-entry observation state is not closed world");
    }
    if (!m_streams->good())
    {
      fail("diagnostic cache stream failure");
    }
  }
  catch (...)
  {
    fail("exception while sealing diagnostic cache");
  }
  return m_error.empty();
}

void DiagnosticWriter::sourceCallback(
    void* user, const c0h2::SourceObservation& observation) noexcept
{
  DiagnosticWriter* writer = static_cast<DiagnosticWriter*>(user);
  if (!writer) return;
  try { writer->writeSource(observation); }
  catch (...) { writer->fail("source observation serialization failed"); }
}

void DiagnosticWriter::scoreCallback(
    void* user, const c0h2::ScoreObservation& observation) noexcept
{
  DiagnosticWriter* writer = static_cast<DiagnosticWriter*>(user);
  if (!writer) return;
  try { writer->writeScore(observation); }
  catch (...) { writer->fail("score observation serialization failed"); }
}

void DiagnosticWriter::eventCallback(
    void* user, const c0rh::EventObservation& observation) noexcept
{
  DiagnosticWriter* writer = static_cast<DiagnosticWriter*>(user);
  if (!writer) return;
  try { writer->writeEvent(observation); }
  catch (...) { writer->fail("event observation serialization failed"); }
}

void DiagnosticWriter::candidateCallback(
    void* user, const c0rh::RawQACandidateObservation& observation) noexcept
{
  DiagnosticWriter* writer = static_cast<DiagnosticWriter*>(user);
  if (!writer) return;
  try { writer->writeCandidate(observation); }
  catch (...) { writer->fail("candidate observation serialization failed"); }
}

void DiagnosticWriter::fillCallback(
    void* user, const c0rh::RawQAFillObservation& observation) noexcept
{
  DiagnosticWriter* writer = static_cast<DiagnosticWriter*>(user);
  if (!writer) return;
  try { writer->writeFill(observation); }
  catch (...) { writer->fail("fill-witness serialization failed"); }
}

void DiagnosticWriter::writeSource(const c0h2::SourceObservation& observation)
{
  std::ostream& out = m_streams->source;
  out << ++m_sequence << '\t'
      << c0h2::toString(observation.kind) << '\t';
  writeOccurrence(out, observation.occurrence);
  out << '\t' << observation.pair.generation << '\t'
      << observation.pair.pair_ordinal << '\t'
      << pointerHex(observation.manager_owner) << '\t'
      << pointerHex(observation.sync_owner) << '\t'
      << pointerHex(observation.io_owner) << '\t'
      << percentEncode(nullable(observation.manager_name)) << '\t'
      << percentEncode(nullable(observation.source_descriptor)) << '\t'
      << observation.source_entry << '\t'
      << observation.read_attempt_ordinal << '\t'
      << c0h2::toString(observation.read_kind) << '\t'
      << c0h2::toString(observation.read_outcome) << '\t'
      << c0h2::toString(observation.phase) << '\t'
      << c0h2::toString(observation.terminal) << '\t'
      << c0h2::toString(observation.transition) << '\t'
      << observation.native_status << '\t'
      << percentEncode(nullable(observation.reason)) << '\n';
  if (!out) throw std::runtime_error("source stream failure");
}

void DiagnosticWriter::writeScore(const c0h2::ScoreObservation& observation)
{
  const std::uint64_t score_sequence = ++m_score_sequence;
  std::ostream& out = m_streams->scores;
  out << score_sequence << '\t';
  writeContext(out, observation.context);
  out << '\t' << percentEncode(nullable(observation.model_reference))
      << '\t' << percentEncode(nullable(observation.scoring_mode))
      << '\t' << percentEncode(nullable(observation.model_sha256))
      << '\t' << percentEncode(nullable(observation.feature_order_sha256))
      << '\t' << percentEncode(nullable(observation.preprocessing_identity))
      << '\t' << percentEncode(nullable(observation.runtime_provider_identity))
      << '\t' << c0h2::toString(observation.status)
      << '\t' << boolean(observation.score_valid)
      << '\t' << (observation.score_valid
                       ? binary32Hex(observation.score) : std::string("NA"))
      << '\t' << observation.failing_feature_index
      << '\t' << observation.feature_count << '\n';

  if (observation.feature_count > 0 &&
      (!observation.ordered_features || !observation.ordered_feature_names))
  {
    throw std::runtime_error("score observation lacks ordered feature payload");
  }
  for (std::size_t index = 0; index < observation.feature_count; ++index)
  {
    m_streams->score_features << score_sequence << '\t' << index << '\t'
        << percentEncode(observation.ordered_feature_names[index]) << '\t'
        << binary32Hex(observation.ordered_features[index]) << '\n';
  }
  if (!out || !m_streams->score_features)
    throw std::runtime_error("score stream failure");
}

void DiagnosticWriter::writeEvent(const c0rh::EventObservation& observation)
{
  const std::uint64_t event_sequence = ++m_sequence;
  std::ostream& out = m_streams->events;
  out << event_sequence << '\t' << c0rh::toString(observation.kind) << '\t'
      << observation.delivery.pair.generation << '\t'
      << observation.delivery.pair.pair_ordinal << '\t'
      << observation.delivery.occurrence_count << '\t';
  writeOccurrence(out, observation.delivery.occurrences[0]);
  out << '\t';
  writeOccurrence(out, observation.delivery.occurrences[1]);
  out << '\t' << observation.delivered_event_ordinal
      << '\t' << boolean(observation.run_valid)
      << '\t' << observation.run_number
      << '\t' << boolean(observation.event_valid)
      << '\t' << observation.event_number
      << '\t' << boolean(observation.centrality_valid)
      << '\t' << binary64Hex(observation.centrality)
      << '\t' << c0rh::toString(observation.outcome)
      << '\t' << observation.native_status
      << '\t' << percentEncode(nullable(observation.outcome_reason))
      << '\t' << observation.encountered_candidate_count
      << '\t' << observation.admitted_candidate_count
      << '\t' << observation.scored_candidate_count
      << '\t' << observation.rawqa_contributor_count
      << '\t' << observation.rawqa_fill_count
      << '\t' << observation.active_trigger_count
      << '\t' << boolean(observation.zero_candidate)
      << '\t' << boolean(observation.zero_fill)
      << '\t' << boolean(observation.finalized)
      << '\t' << percentEncode(nullable(observation.module_name)) << '\n';
  if (observation.active_trigger_count > 0 && !observation.active_triggers)
    throw std::runtime_error("event observation lacks active-trigger payload");
  for (std::size_t index = 0; index < observation.active_trigger_count; ++index)
    m_streams->event_triggers << event_sequence << '\t' << index << '\t'
        << percentEncode(observation.active_triggers[index]) << '\n';
  if (!out || !m_streams->event_triggers)
    throw std::runtime_error("event stream failure");
}

void DiagnosticWriter::writeCandidate(
    const c0rh::RawQACandidateObservation& observation)
{
  const std::uint64_t candidate_sequence = ++m_candidate_sequence;
  std::ostream& out = m_streams->candidates;
  out << candidate_sequence << '\t';
  writeContext(out, observation.context);
  out << '\t' << observation.photon_pt_slice
      << '\t' << observation.centrality_slice
      << '\t' << boolean(observation.canonical_view)
      << '\t' << percentEncode(nullable(observation.view_suffix))
      << '\t' << c0rh::toString(observation.tight_tag)
      << '\t' << boolean(observation.preselection_pass)
      << '\t' << boolean(observation.direct_candidate_admitted)
      << '\t' << observation.valid_mask
      << '\t' << observation.value_count
      << '\t' << observation.active_trigger_count << '\n';

  if (observation.value_count > 0 &&
      (!observation.values || !observation.value_bits ||
       !observation.variable_names))
  {
    throw std::runtime_error("candidate observation lacks value payload");
  }
  for (std::size_t index = 0; index < observation.value_count; ++index)
  {
    if (binary64Hex(observation.values[index]) !=
        unsignedHex(observation.value_bits[index], 16))
    {
      throw std::runtime_error("candidate binary64 value/bits disagreement");
    }
    m_streams->candidate_values << candidate_sequence << '\t' << index << '\t'
        << percentEncode(observation.variable_names[index]) << '\t'
        << unsignedHex(observation.value_bits[index], 16) << '\t'
        << boolean((observation.valid_mask & (1U << index)) != 0) << '\n';
  }
  if (observation.active_trigger_count > 0 && !observation.active_triggers)
  {
    throw std::runtime_error("candidate observation lacks trigger payload");
  }
  for (std::size_t index = 0; index < observation.active_trigger_count; ++index)
  {
    m_streams->candidate_triggers << candidate_sequence << '\t' << index << '\t'
        << percentEncode(observation.active_triggers[index]) << '\n';
  }
  if (!out || !m_streams->candidate_values || !m_streams->candidate_triggers)
    throw std::runtime_error("candidate stream failure");
}

void DiagnosticWriter::writeFill(const c0rh::RawQAFillObservation& observation)
{
  std::ostream& out = m_streams->fills;
  out << ++m_sequence << '\t';
  writeContext(out, observation.context);
  out << '\t' << percentEncode(nullable(observation.variable_name))
      << '\t' << percentEncode(nullable(observation.trigger_name))
      << '\t' << percentEncode(nullable(observation.tag_name))
      << '\t' << percentEncode(nullable(observation.view_suffix))
      << '\t' << percentEncode(nullable(observation.directory_name))
      << '\t' << percentEncode(nullable(observation.object_name))
      << '\t' << percentEncode(nullable(observation.object_class))
      << '\t' << percentEncode(nullable(observation.object_contract_key))
      << '\t' << observation.photon_pt_slice
      << '\t' << observation.centrality_slice
      << '\t' << boolean(observation.canonical_view)
      << '\t' << c0rh::toString(observation.tight_tag)
      << '\t' << unsignedHex(observation.value_bits, 16)
      << '\t' << unsignedHex(observation.weight_bits, 16)
      << '\t' << observation.fill_ordinal
      << '\t' << boolean(observation.value_valid)
      << '\t' << boolean(observation.filled)
      << '\t' << boolean(observation.sumw2_enabled) << '\n';
  if (binary64Hex(observation.value) != unsignedHex(observation.value_bits, 16) ||
      binary64Hex(observation.weight) != unsignedHex(observation.weight_bits, 16))
  {
    throw std::runtime_error("fill witness binary64 value/bits disagreement");
  }
  if (!out) throw std::runtime_error("fill-witness stream failure");
}

}  // namespace c0r
}  // namespace the106
