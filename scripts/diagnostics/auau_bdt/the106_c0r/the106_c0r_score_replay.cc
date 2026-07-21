// Cache-only same-runtime score closure for AU-AU-PHOTON-C0-RAW-QA.
#include "the106_c0r_diagnostic_writer.h"

#include <TMVA/RBDT.hxx>
#include <TROOT.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

typedef std::vector<std::string> Row;

struct Table
{
  std::map<std::string, std::size_t> columns;
  std::vector<Row> rows;
};

Row split(const std::string& line)
{
  Row result;
  std::size_t begin = 0;
  for (std::size_t index = 0; index <= line.size(); ++index)
  {
    if (index == line.size() || line[index] == '\t')
    {
      result.push_back(line.substr(begin, index - begin));
      begin = index + 1;
    }
  }
  return result;
}

Table readTable(const std::string& path, const Row& expected_header)
{
  std::ifstream input(path.c_str());
  if (!input) throw std::runtime_error("cannot open required table " + path);
  std::string line;
  if (!std::getline(input, line)) throw std::runtime_error("empty table " + path);
  if (!line.empty() && line.back() == '\r') line.pop_back();
  const Row header = split(line);
  if (header != expected_header)
    throw std::runtime_error("table schema mismatch in " + path);
  Table table;
  for (std::size_t index = 0; index < header.size(); ++index)
  {
    if (header[index].empty() ||
        !table.columns.insert(std::make_pair(header[index], index)).second)
      throw std::runtime_error("duplicate or empty table column in " + path);
  }
  while (std::getline(input, line))
  {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty())
      throw std::runtime_error("blank row in " + path);
    const Row row = split(line);
    if (row.size() != header.size())
      throw std::runtime_error("table width mismatch in " + path);
    table.rows.push_back(row);
  }
  if (!input.eof()) throw std::runtime_error("read failure in " + path);
  return table;
}

const std::string& field(const Table& table, const Row& row, const char* name)
{
  const std::map<std::string, std::size_t>::const_iterator found =
      table.columns.find(name);
  if (found == table.columns.end())
    throw std::runtime_error(std::string("missing required column ") + name);
  return row.at(found->second);
}

std::uint64_t unsignedValue(const std::string& text, const char* label)
{
  if (text.empty() ||
      !std::all_of(text.begin(), text.end(),
                   [](char ch) { return ch >= '0' && ch <= '9'; }))
    throw std::runtime_error(std::string("invalid unsigned value for ") + label);
  std::size_t used = 0;
  const unsigned long long result = std::stoull(text, &used, 10);
  if (used != text.size())
    throw std::runtime_error(std::string("invalid unsigned value for ") + label);
  return static_cast<std::uint64_t>(result);
}

bool booleanValue(const std::string& text, const char* label)
{
  if (text == "0") return false;
  if (text == "1") return true;
  throw std::runtime_error(std::string("invalid boolean value for ") + label);
}

bool isLowerSha256(const std::string& value)
{
  if (value.size() != 64) return false;
  for (std::string::const_iterator iter = value.begin(); iter != value.end(); ++iter)
    if (!((*iter >= '0' && *iter <= '9') || (*iter >= 'a' && *iter <= 'f')))
      return false;
  return true;
}

// Small self-contained SHA-256 implementation.  It keeps model-byte authority
// local to this bounded validator and avoids accepting a caller-authored hash.
class Sha256
{
 public:
  Sha256()
    : m_state{{0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U, 0xa54ff53aU,
               0x510e527fU, 0x9b05688cU, 0x1f83d9abU, 0x5be0cd19U}}
    , m_buffer{{}}
    , m_buffer_size(0)
    , m_total_size(0)
  {
  }

  void update(const unsigned char* data, std::size_t size)
  {
    if (!data && size != 0) throw std::runtime_error("null SHA-256 input");
    m_total_size += static_cast<std::uint64_t>(size);
    while (size > 0)
    {
      const std::size_t copied = std::min(size, m_buffer.size() - m_buffer_size);
      std::memcpy(m_buffer.data() + m_buffer_size, data, copied);
      m_buffer_size += copied;
      data += copied;
      size -= copied;
      if (m_buffer_size == m_buffer.size())
      {
        transform(m_buffer.data());
        m_buffer_size = 0;
      }
    }
  }

  void update(const std::string& value)
  {
    update(reinterpret_cast<const unsigned char*>(value.data()), value.size());
  }

  std::string final()
  {
    const std::uint64_t bit_count = m_total_size * 8U;
    unsigned char one = 0x80U;
    updatePadding(&one, 1);
    const unsigned char zero = 0;
    while (m_buffer_size != 56) updatePadding(&zero, 1);
    unsigned char length[8];
    for (std::size_t index = 0; index < 8; ++index)
      length[7 - index] = static_cast<unsigned char>((bit_count >> (8U * index)) & 0xffU);
    updatePadding(length, sizeof(length));
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (std::size_t index = 0; index < m_state.size(); ++index)
      out << std::setw(8) << m_state[index];
    return out.str();
  }

 private:
  static std::uint32_t rotate(std::uint32_t value, unsigned int count)
  {
    return (value >> count) | (value << (32U - count));
  }

  void updatePadding(const unsigned char* data, std::size_t size)
  {
    while (size > 0)
    {
      const std::size_t copied = std::min(size, m_buffer.size() - m_buffer_size);
      std::memcpy(m_buffer.data() + m_buffer_size, data, copied);
      m_buffer_size += copied;
      data += copied;
      size -= copied;
      if (m_buffer_size == m_buffer.size())
      {
        transform(m_buffer.data());
        m_buffer_size = 0;
      }
    }
  }

  void transform(const unsigned char* block)
  {
    static const std::uint32_t constants[64] = {
        0x428a2f98U,0x71374491U,0xb5c0fbcfU,0xe9b5dba5U,0x3956c25bU,0x59f111f1U,0x923f82a4U,0xab1c5ed5U,
        0xd807aa98U,0x12835b01U,0x243185beU,0x550c7dc3U,0x72be5d74U,0x80deb1feU,0x9bdc06a7U,0xc19bf174U,
        0xe49b69c1U,0xefbe4786U,0x0fc19dc6U,0x240ca1ccU,0x2de92c6fU,0x4a7484aaU,0x5cb0a9dcU,0x76f988daU,
        0x983e5152U,0xa831c66dU,0xb00327c8U,0xbf597fc7U,0xc6e00bf3U,0xd5a79147U,0x06ca6351U,0x14292967U,
        0x27b70a85U,0x2e1b2138U,0x4d2c6dfcU,0x53380d13U,0x650a7354U,0x766a0abbU,0x81c2c92eU,0x92722c85U,
        0xa2bfe8a1U,0xa81a664bU,0xc24b8b70U,0xc76c51a3U,0xd192e819U,0xd6990624U,0xf40e3585U,0x106aa070U,
        0x19a4c116U,0x1e376c08U,0x2748774cU,0x34b0bcb5U,0x391c0cb3U,0x4ed8aa4aU,0x5b9cca4fU,0x682e6ff3U,
        0x748f82eeU,0x78a5636fU,0x84c87814U,0x8cc70208U,0x90befffaU,0xa4506cebU,0xbef9a3f7U,0xc67178f2U};
    std::uint32_t words[64];
    for (std::size_t index = 0; index < 16; ++index)
    {
      const std::size_t offset = 4 * index;
      words[index] = (static_cast<std::uint32_t>(block[offset]) << 24U) |
                     (static_cast<std::uint32_t>(block[offset + 1]) << 16U) |
                     (static_cast<std::uint32_t>(block[offset + 2]) << 8U) |
                     static_cast<std::uint32_t>(block[offset + 3]);
    }
    for (std::size_t index = 16; index < 64; ++index)
    {
      const std::uint32_t s0 = rotate(words[index - 15], 7) ^
                               rotate(words[index - 15], 18) ^
                               (words[index - 15] >> 3U);
      const std::uint32_t s1 = rotate(words[index - 2], 17) ^
                               rotate(words[index - 2], 19) ^
                               (words[index - 2] >> 10U);
      words[index] = words[index - 16] + s0 + words[index - 7] + s1;
    }
    std::uint32_t a = m_state[0], b = m_state[1], c = m_state[2], d = m_state[3];
    std::uint32_t e = m_state[4], f = m_state[5], g = m_state[6], h = m_state[7];
    for (std::size_t index = 0; index < 64; ++index)
    {
      const std::uint32_t s1 = rotate(e, 6) ^ rotate(e, 11) ^ rotate(e, 25);
      const std::uint32_t choice = (e & f) ^ ((~e) & g);
      const std::uint32_t t1 = h + s1 + choice + constants[index] + words[index];
      const std::uint32_t s0 = rotate(a, 2) ^ rotate(a, 13) ^ rotate(a, 22);
      const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t t2 = s0 + majority;
      h = g; g = f; f = e; e = d + t1;
      d = c; c = b; b = a; a = t1 + t2;
    }
    m_state[0] += a; m_state[1] += b; m_state[2] += c; m_state[3] += d;
    m_state[4] += e; m_state[5] += f; m_state[6] += g; m_state[7] += h;
  }

  std::array<std::uint32_t, 8> m_state;
  std::array<unsigned char, 64> m_buffer;
  std::size_t m_buffer_size;
  std::uint64_t m_total_size;
};

std::string sha256File(const std::string& path)
{
  std::ifstream input(path.c_str(), std::ios::in | std::ios::binary);
  if (!input) throw std::runtime_error("cannot open sealed model file " + path);
  Sha256 hash;
  std::array<char, 65536> buffer;
  while (input)
  {
    input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const std::streamsize count = input.gcount();
    if (count > 0)
      hash.update(reinterpret_cast<const unsigned char*>(buffer.data()),
                  static_cast<std::size_t>(count));
  }
  if (!input.eof()) throw std::runtime_error("read failure in sealed model file " + path);
  return hash.final();
}

std::string sha256FeatureOrder(const std::vector<std::string>& names)
{
  Sha256 hash;
  for (std::vector<std::string>::const_iterator iter = names.begin();
       iter != names.end(); ++iter)
  {
    if (iter->empty() ||
        !std::all_of(iter->begin(), iter->end(), [](char ch) {
          return (ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') ||
                 (ch >= '0' && ch <= '9') || ch == '_' || ch == '.' || ch == '-';
        }))
      throw std::runtime_error("unsafe or empty feature name in ordered feature contract");
    hash.update(*iter);
    hash.update("\n");
  }
  return hash.final();
}

std::string canonicalRuntimeIdentity()
{
  if (!gROOT || !gROOT->GetVersion())
    throw std::runtime_error("ROOT runtime identity is unavailable");
  return std::string("ROOT/") + gROOT->GetVersion() +
         "|TMVA::Experimental::RBDT|model_key=myBDT";
}

std::string jsonEscape(const std::string& value)
{
  std::ostringstream out;
  for (std::string::const_iterator iter = value.begin(); iter != value.end(); ++iter)
  {
    const unsigned char ch = static_cast<unsigned char>(*iter);
    switch (ch)
    {
      case '"': out << "\\\""; break;
      case '\\': out << "\\\\"; break;
      case '\b': out << "\\b"; break;
      case '\f': out << "\\f"; break;
      case '\n': out << "\\n"; break;
      case '\r': out << "\\r"; break;
      case '\t': out << "\\t"; break;
      default:
        if (ch < 0x20U)
          out << "\\u" << std::hex << std::setfill('0') << std::setw(4)
              << static_cast<unsigned int>(ch) << std::dec;
        else
          out << static_cast<char>(ch);
    }
  }
  return out.str();
}

void writeReport(const std::string& path, const std::string& body)
{
  std::ofstream report(path.c_str(), std::ios::out | std::ios::trunc);
  if (!report) throw std::runtime_error("cannot create deterministic score report " + path);
  report << body;
  if (!report) throw std::runtime_error("cannot write deterministic score report " + path);
}

std::string joinKey(const Table& table, const Row& row,
                    const std::vector<const char*>& columns)
{
  std::string result;
  for (std::size_t index = 0; index < columns.size(); ++index)
  {
    if (index != 0) result.push_back('\x1f');
    const std::string& value = field(table, row, columns[index]);
    if (value.find('\x1f') != std::string::npos)
      throw std::runtime_error("unit separator in identity field");
    result += value;
  }
  return result;
}

std::string eventKey(const Table& table, const Row& row)
{
  static const std::vector<const char*> columns = {
      "pair_generation", "pair_ordinal", "delivered_event_ordinal",
      "run_valid", "run_number", "event_valid", "event_number"};
  return joinKey(table, row, columns);
}

std::string candidateKey(const Table& table, const Row& row, bool score_table)
{
  const std::vector<const char*> columns = {
      "pair_generation", "pair_ordinal", "delivered_event_ordinal",
      "run_valid", "run_number", "event_valid", "event_number",
      "container_key", "cluster_id", "producer_encounter_ordinal",
      score_table ? "canonical_view" : "canonical_view_context",
      "view_key", "module_name"};
  return joinKey(table, row, columns);
}

struct Feature
{
  std::string name;
  std::string bits;
};

struct Summary
{
  std::uint64_t observed_rows;
  std::uint64_t compute_rows;
  std::uint64_t valid_rows;
  std::uint64_t exact_rows;
  std::uint64_t noncompute_rows;
  std::uint64_t event_rows;
  std::uint64_t candidate_rows;
  std::string model_reference;
  std::string model_sha256;
  std::string feature_order_sha256;
  std::string preprocessing_identity;
  std::string runtime_provider_identity;
  std::string scoring_mode;
  std::vector<std::string> feature_names;
};

Summary closeScores(const std::string& cache)
{
  const Row metadata_header = {"key", "value"};
  const Row scores_header = {
      "sequence", "pair_generation", "pair_ordinal", "delivered_event_ordinal",
      "run_valid", "run_number", "event_valid", "event_number", "container_key",
      "cluster_id", "producer_encounter_ordinal", "canonical_view", "view_key",
      "module_name", "model_reference", "scoring_mode", "model_sha256",
      "feature_order_sha256", "preprocessing_identity", "runtime_provider_identity",
      "status", "score_valid", "score_bits", "failing_feature_index", "feature_count"};
  const Row features_header = {
      "score_sequence", "feature_index", "feature_name", "value_bits"};
  const Row events_header = {
      "sequence", "kind", "pair_generation", "pair_ordinal", "occurrence_count",
      "occurrence0_generation", "occurrence0_manager_ordinal", "occurrence0_ordinal",
      "occurrence1_generation", "occurrence1_manager_ordinal", "occurrence1_ordinal",
      "delivered_event_ordinal", "run_valid", "run_number", "event_valid", "event_number",
      "centrality_valid", "centrality_bits", "outcome", "native_status", "outcome_reason",
      "encountered_candidate_count", "admitted_candidate_count", "scored_candidate_count",
      "rawqa_contributor_count", "rawqa_fill_count", "active_trigger_count",
      "zero_candidate", "zero_fill", "finalized", "module_name"};
  const Row candidates_header = {
      "sequence", "pair_generation", "pair_ordinal", "delivered_event_ordinal",
      "run_valid", "run_number", "event_valid", "event_number", "container_key",
      "cluster_id", "producer_encounter_ordinal", "canonical_view_context", "view_key",
      "module_name", "photon_pt_slice", "centrality_slice", "canonical_view", "view_suffix",
      "tight_tag", "preselection_pass", "direct_candidate_admitted", "valid_mask",
      "value_count", "active_trigger_count"};

  const Table metadata = readTable(cache + "/metadata.tsv", metadata_header);
  const Table scores = readTable(cache + "/score_observations.tsv", scores_header);
  const Table features = readTable(cache + "/score_features.tsv", features_header);
  const Table events = readTable(cache + "/event_observations.tsv", events_header);
  const Table candidates = readTable(cache + "/rawqa_candidates.tsv", candidates_header);

  static const char* required_metadata[] = {
      "format", "canary_id", "scientific_contract_sha256", "source_manifest_sha256",
      "configuration_sha256", "reconstruction_identity", "cdb_identity",
      "model_reference", "model_sha256", "feature_order_sha256",
      "preprocessing_identity", "runtime_provider_identity", "source_range"};
  std::map<std::string, std::string> metadata_values;
  for (std::vector<Row>::const_iterator row = metadata.rows.begin();
       row != metadata.rows.end(); ++row)
  {
    const std::string key = field(metadata, *row, "key");
    const std::string value = the106::c0r::percentDecode(field(metadata, *row, "value"));
    if (!metadata_values.insert(std::make_pair(key, value)).second)
      throw std::runtime_error("duplicate metadata key " + key);
  }
  if (metadata_values.size() !=
      sizeof(required_metadata) / sizeof(required_metadata[0]))
    throw std::runtime_error("unexpected metadata key set");
  for (std::size_t index = 0;
       index < sizeof(required_metadata) / sizeof(required_metadata[0]); ++index)
    if (metadata_values.find(required_metadata[index]) == metadata_values.end() ||
        metadata_values[required_metadata[index]].empty())
      throw std::runtime_error(std::string("missing metadata authority ") +
                               required_metadata[index]);
  if (metadata_values["format"] != "THE106_C0R_CACHE_V1")
    throw std::runtime_error("unsupported cache format");
  if (!isLowerSha256(metadata_values["model_sha256"]) ||
      !isLowerSha256(metadata_values["feature_order_sha256"]))
    throw std::runtime_error("malformed model or feature-order SHA-256 authority");

  const std::string runtime_identity = canonicalRuntimeIdentity();
  if (metadata_values["runtime_provider_identity"] != runtime_identity)
    throw std::runtime_error("runtime/provider identity does not match executing ROOT/TMVA runtime");
  const std::string actual_model_sha = sha256File(metadata_values["model_reference"]);
  if (actual_model_sha != metadata_values["model_sha256"])
    throw std::runtime_error("sealed model bytes SHA-256 mismatch");

  std::map<std::uint64_t, std::map<std::size_t, Feature> > by_score;
  std::uint64_t expected_feature_row = 0;
  for (std::vector<Row>::const_iterator row = features.rows.begin();
       row != features.rows.end(); ++row)
  {
    ++expected_feature_row;
    const std::uint64_t sequence = unsignedValue(
        field(features, *row, "score_sequence"), "score_sequence");
    const std::size_t index = static_cast<std::size_t>(unsignedValue(
        field(features, *row, "feature_index"), "feature_index"));
    Feature feature;
    feature.name = the106::c0r::percentDecode(field(features, *row, "feature_name"));
    feature.bits = field(features, *row, "value_bits");
    float value = 0.0F;
    if (!the106::c0r::binary32FromHex(feature.bits, value) ||
        the106::c0r::binary32Hex(value) != feature.bits)
      throw std::runtime_error("invalid or noncanonical feature bits");
    if (!by_score[sequence].insert(std::make_pair(index, feature)).second)
      throw std::runtime_error("duplicate feature index in score observation");
    (void) expected_feature_row;
  }

  std::map<std::string, std::uint64_t> event_scored_expected;
  std::map<std::string, std::uint64_t> event_scored_observed;
  std::set<std::string> finalized_events;
  for (std::vector<Row>::const_iterator row = events.rows.begin();
       row != events.rows.end(); ++row)
  {
    if (field(events, *row, "kind") != "EventObservationKind::finalization") continue;
    if (!booleanValue(field(events, *row, "finalized"), "finalized"))
      throw std::runtime_error("event finalization row is not finalized");
    const std::string key = eventKey(events, *row);
    if (!finalized_events.insert(key).second)
      throw std::runtime_error("duplicate event finalization identity");
    event_scored_expected[key] = unsignedValue(
        field(events, *row, "scored_candidate_count"), "scored_candidate_count");
  }

  std::set<std::string> admitted_candidates;
  for (std::vector<Row>::const_iterator row = candidates.rows.begin();
       row != candidates.rows.end(); ++row)
  {
    if (!booleanValue(field(candidates, *row, "direct_candidate_admitted"),
                      "direct_candidate_admitted"))
      throw std::runtime_error("raw-QA cache contains a nonadmitted candidate row");
    if (!booleanValue(field(candidates, *row, "canonical_view_context"),
                      "canonical_view_context") ||
        !booleanValue(field(candidates, *row, "canonical_view"), "canonical_view"))
      throw std::runtime_error("raw-QA cache contains a noncanonical candidate row");
    const std::string key = candidateKey(candidates, *row, false);
    if (!admitted_candidates.insert(key).second)
      throw std::runtime_error("duplicate admitted raw-QA candidate identity");
    if (finalized_events.find(eventKey(candidates, *row)) == finalized_events.end())
      throw std::runtime_error("admitted raw-QA candidate lacks event finalization");
  }

  Summary summary = {};
  summary.model_reference = metadata_values["model_reference"];
  summary.model_sha256 = actual_model_sha;
  summary.feature_order_sha256 = metadata_values["feature_order_sha256"];
  summary.preprocessing_identity = metadata_values["preprocessing_identity"];
  summary.runtime_provider_identity = runtime_identity;
  summary.event_rows = finalized_events.size();
  summary.candidate_rows = admitted_candidates.size();

  std::set<std::uint64_t> sequences;
  std::set<std::string> score_candidate_keys;
  std::set<std::string> scored_canonical_candidates;
  std::vector<std::string> frozen_feature_names;
  std::string scoring_mode;
  std::unique_ptr<TMVA::Experimental::RBDT> model;
  for (std::size_t row_index = 0; row_index < scores.rows.size(); ++row_index)
  {
    const Row& row = scores.rows[row_index];
    const std::uint64_t sequence = unsignedValue(field(scores, row, "sequence"), "sequence");
    if (sequence != row_index + 1)
      throw std::runtime_error("score sequence is not contiguous writer order");
    if (!sequences.insert(sequence).second)
      throw std::runtime_error("duplicate score sequence");
    ++summary.observed_rows;

    const std::string row_model = the106::c0r::percentDecode(
        field(scores, row, "model_reference"));
    const std::string row_model_sha = the106::c0r::percentDecode(
        field(scores, row, "model_sha256"));
    const std::string row_feature_sha = the106::c0r::percentDecode(
        field(scores, row, "feature_order_sha256"));
    const std::string row_preprocessing = the106::c0r::percentDecode(
        field(scores, row, "preprocessing_identity"));
    const std::string row_runtime = the106::c0r::percentDecode(
        field(scores, row, "runtime_provider_identity"));
    const std::string row_mode = the106::c0r::percentDecode(
        field(scores, row, "scoring_mode"));
    if (row_model != summary.model_reference || row_model_sha != summary.model_sha256 ||
        row_feature_sha != summary.feature_order_sha256 ||
        row_preprocessing != summary.preprocessing_identity ||
        row_runtime != summary.runtime_provider_identity)
      throw std::runtime_error("score-row provenance differs from sealed metadata authority");
    if (row_mode.empty()) throw std::runtime_error("empty score mode");
    if (scoring_mode.empty()) scoring_mode = row_mode;
    if (scoring_mode != row_mode)
      throw std::runtime_error("multiple score modes in single-model C0-R cache");

    const std::string score_candidate = candidateKey(scores, row, true);
    if (!score_candidate_keys.insert(score_candidate).second)
      throw std::runtime_error("duplicate score observation for candidate/view identity");
    const std::string score_event = eventKey(scores, row);
    if (finalized_events.find(score_event) == finalized_events.end())
      throw std::runtime_error("score observation lacks event finalization");

    const std::size_t count = static_cast<std::size_t>(unsignedValue(
        field(scores, row, "feature_count"), "feature_count"));
    std::vector<float> input;
    std::vector<std::string> names;
    input.reserve(count);
    names.reserve(count);
    const std::map<std::uint64_t, std::map<std::size_t, Feature> >::const_iterator group =
        by_score.find(sequence);
    if (count > 0 && group == by_score.end())
      throw std::runtime_error("missing ordered feature group");
    const std::map<std::size_t, Feature> empty_features;
    const std::map<std::size_t, Feature>& feature_map =
        group == by_score.end() ? empty_features : group->second;
    if (feature_map.size() != count)
      throw std::runtime_error("ordered feature cardinality mismatch");
    for (std::size_t index = 0; index < count; ++index)
    {
      const std::map<std::size_t, Feature>::const_iterator feature =
          feature_map.find(index);
      if (feature == feature_map.end())
        throw std::runtime_error("missing ordered feature index");
      float value = 0.0F;
      if (!the106::c0r::binary32FromHex(feature->second.bits, value))
        throw std::runtime_error("invalid ordered feature bits");
      input.push_back(value);
      names.push_back(feature->second.name);
    }

    const std::string status = field(scores, row, "status");
    const bool score_valid = booleanValue(field(scores, row, "score_valid"), "score_valid");
    const std::size_t failing_index = static_cast<std::size_t>(unsignedValue(
        field(scores, row, "failing_feature_index"), "failing_feature_index"));
    if (status == "ScoreStatus::feature_nonfinite")
    {
      if (score_valid || field(scores, row, "score_bits") != "NA" ||
          failing_index != count)
        throw std::runtime_error("feature-nonfinite score-state contract mismatch");
      ++summary.noncompute_rows;
      continue;
    }
    if (status != "ScoreStatus::valid" && status != "ScoreStatus::output_empty" &&
        status != "ScoreStatus::output_nonfinite")
      throw std::runtime_error("unsupported score status");
    if (failing_index != 0)
      throw std::runtime_error("computed score row has nonzero failing feature index");
    if (frozen_feature_names.empty()) frozen_feature_names = names;
    if (names != frozen_feature_names)
      throw std::runtime_error("ordered feature names differ between Compute rows");
    if (sha256FeatureOrder(names) != summary.feature_order_sha256)
      throw std::runtime_error("observed ordered feature names do not match feature-order SHA-256");

    if (!model)
      model.reset(new TMVA::Experimental::RBDT("myBDT", summary.model_reference));
    const std::vector<float> output = model->Compute(input);
    ++summary.compute_rows;
    if (status == "ScoreStatus::output_empty")
    {
      if (score_valid || field(scores, row, "score_bits") != "NA" || !output.empty())
        throw std::runtime_error("output-empty score-state changed");
      ++summary.exact_rows;
    }
    else
    {
      if (!score_valid || output.empty())
        throw std::runtime_error("computed score lacks direct or replay output");
      float direct_score = 0.0F;
      if (!the106::c0r::binary32FromHex(field(scores, row, "score_bits"), direct_score) ||
          the106::c0r::binary32Hex(direct_score) != field(scores, row, "score_bits"))
        throw std::runtime_error("invalid or noncanonical direct score bits");
      if (the106::c0r::binary32Hex(output[0]) != field(scores, row, "score_bits"))
        throw std::runtime_error("same-runtime score bits differ");
      if ((status == "ScoreStatus::valid") != std::isfinite(output[0]))
        throw std::runtime_error("same-runtime score finite state differs");
      if (status == "ScoreStatus::valid") ++summary.valid_rows;
      ++summary.exact_rows;
    }

    if (booleanValue(field(scores, row, "canonical_view"), "canonical_view"))
    {
      scored_canonical_candidates.insert(score_candidate);
      event_scored_observed[score_event] =
          static_cast<std::uint64_t>(std::count_if(
              scored_canonical_candidates.begin(), scored_canonical_candidates.end(),
              [&score_event](const std::string& key) {
                return key.compare(0, score_event.size(), score_event) == 0 &&
                       key.size() > score_event.size() && key[score_event.size()] == '\x1f';
              }));
    }
  }

  for (std::map<std::uint64_t, std::map<std::size_t, Feature> >::const_iterator row =
           by_score.begin(); row != by_score.end(); ++row)
    if (sequences.find(row->first) == sequences.end())
      throw std::runtime_error("orphan score feature group");
  if (summary.compute_rows == 0 || summary.valid_rows == 0)
    throw std::runtime_error("vacuous score closure: no valid Compute observation");
  for (std::set<std::string>::const_iterator candidate = admitted_candidates.begin();
       candidate != admitted_candidates.end(); ++candidate)
    if (score_candidate_keys.find(*candidate) == score_candidate_keys.end())
      throw std::runtime_error("admitted raw-QA candidate lacks a score observation");
  for (std::map<std::string, std::uint64_t>::const_iterator event =
           event_scored_expected.begin(); event != event_scored_expected.end(); ++event)
  {
    const std::uint64_t observed = event_scored_observed[event->first];
    if (observed != event->second)
      throw std::runtime_error("event scored-candidate count does not close");
  }

  summary.scoring_mode = scoring_mode;
  summary.feature_names = frozen_feature_names;
  return summary;
}

std::string passJson(const Summary& summary)
{
  std::ostringstream out;
  out << "{\n"
      << "  \"result\": \"PASS\",\n"
      << "  \"comparison\": \"EXACT_SAME_RUNTIME_RBDT_COMPUTE\",\n"
      << "  \"model_key\": \"myBDT\",\n"
      << "  \"model_reference\": \"" << jsonEscape(summary.model_reference) << "\",\n"
      << "  \"model_sha256\": \"" << summary.model_sha256 << "\",\n"
      << "  \"feature_order_sha256\": \"" << summary.feature_order_sha256 << "\",\n"
      << "  \"preprocessing_identity\": \""
      << jsonEscape(summary.preprocessing_identity) << "\",\n"
      << "  \"runtime_provider_identity\": \""
      << jsonEscape(summary.runtime_provider_identity) << "\",\n"
      << "  \"scoring_mode\": \"" << jsonEscape(summary.scoring_mode) << "\",\n"
      << "  \"event_finalizations\": " << summary.event_rows << ",\n"
      << "  \"admitted_rawqa_candidates\": " << summary.candidate_rows << ",\n"
      << "  \"observed_score_rows\": " << summary.observed_rows << ",\n"
      << "  \"computed_score_rows\": " << summary.compute_rows << ",\n"
      << "  \"valid_score_rows\": " << summary.valid_rows << ",\n"
      << "  \"exact_score_rows\": " << summary.exact_rows << ",\n"
      << "  \"noncompute_rows\": " << summary.noncompute_rows << ",\n"
      << "  \"ordered_feature_names\": [";
  for (std::size_t index = 0; index < summary.feature_names.size(); ++index)
  {
    if (index != 0) out << ", ";
    out << "\"" << jsonEscape(summary.feature_names[index]) << "\"";
  }
  out << "]\n}\n";
  return out.str();
}

std::string failJson(const std::string& error)
{
  return std::string("{\n  \"result\": \"FAIL\",\n") +
         "  \"comparison\": \"EXACT_SAME_RUNTIME_RBDT_COMPUTE\",\n" +
         "  \"error\": \"" + jsonEscape(error) + "\"\n}\n";
}

}  // namespace

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    std::cerr << "usage: " << argv[0] << " CACHE_DIR REPORT.json\n";
    return 2;
  }
  try
  {
    const Summary summary = closeScores(argv[1]);
    writeReport(argv[2], passJson(summary));
    return 0;
  }
  catch (const std::exception& error)
  {
    try { writeReport(argv[2], failJson(error.what())); }
    catch (...) {}
    std::cerr << "THE106 C0-R score closure fail-closed: " << error.what() << '\n';
    return 1;
  }
}
