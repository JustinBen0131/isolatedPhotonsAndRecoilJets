#include "the106_c0r_diagnostic_writer.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TKey.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
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

struct Table
{
  std::vector<std::string> header;
  std::map<std::string, std::size_t> columns;
  std::vector<std::vector<std::string> > rows;
};

std::vector<std::string> split(const std::string& line, char separator)
{
  std::vector<std::string> fields;
  std::size_t begin = 0;
  for (std::size_t index = 0; index <= line.size(); ++index)
  {
    if (index == line.size() || line[index] == separator)
    {
      fields.push_back(line.substr(begin, index - begin));
      begin = index + 1;
    }
  }
  return fields;
}

Table readTable(const std::string& path)
{
  std::ifstream input(path.c_str());
  if (!input) throw std::runtime_error("cannot open table " + path);
  Table table;
  std::string line;
  if (!std::getline(input, line)) throw std::runtime_error("empty table " + path);
  if (!line.empty() && line[line.size() - 1] == '\r') line.resize(line.size() - 1);
  table.header = split(line, '\t');
  for (std::size_t index = 0; index < table.header.size(); ++index)
  {
    if (!table.columns.insert(std::make_pair(table.header[index], index)).second)
      throw std::runtime_error("duplicate table column " + table.header[index]);
  }
  while (std::getline(input, line))
  {
    if (!line.empty() && line[line.size() - 1] == '\r') line.resize(line.size() - 1);
    if (line.empty()) continue;
    std::vector<std::string> row = split(line, '\t');
    if (row.size() != table.header.size())
      throw std::runtime_error("inconsistent row width in " + path);
    table.rows.push_back(row);
  }
  return table;
}

const std::string& field(const Table& table,
                         const std::vector<std::string>& row,
                         const std::string& name)
{
  const std::map<std::string, std::size_t>::const_iterator found =
      table.columns.find(name);
  if (found == table.columns.end())
    throw std::runtime_error("missing table column " + name);
  return row.at(found->second);
}

std::uint64_t unsignedValue(const std::string& text)
{
  std::size_t used = 0;
  const unsigned long long value = std::stoull(text, &used, 10);
  if (used != text.size()) throw std::runtime_error("invalid unsigned integer " + text);
  return static_cast<std::uint64_t>(value);
}

int integerValue(const std::string& text)
{
  std::size_t used = 0;
  const long value = std::stol(text, &used, 10);
  if (used != text.size()) throw std::runtime_error("invalid integer " + text);
  return static_cast<int>(value);
}

bool booleanValue(const std::string& text)
{
  if (text == "1") return true;
  if (text == "0") return false;
  throw std::runtime_error("invalid boolean " + text);
}

std::string jsonEscape(const std::string& value)
{
  std::ostringstream out;
  for (std::string::const_iterator iter = value.begin(); iter != value.end(); ++iter)
  {
    const unsigned char ch = static_cast<unsigned char>(*iter);
    if (ch == '\\' || ch == '"') out << '\\' << static_cast<char>(ch);
    else if (ch == '\n') out << "\\n";
    else if (ch == '\r') out << "\\r";
    else if (ch == '\t') out << "\\t";
    else if (ch < 0x20)
      out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
          << static_cast<unsigned int>(ch) << std::dec;
    else out << static_cast<char>(ch);
  }
  return out.str();
}

std::string contextKey(const Table& table, const std::vector<std::string>& row)
{
  const char* names[] = {
      "pair_generation", "pair_ordinal", "delivered_event_ordinal",
      "run_valid", "run_number", "event_valid", "event_number",
      "container_key", "cluster_id", "producer_encounter_ordinal",
      "canonical_view_context", "view_key", "module_name"};
  std::ostringstream out;
  for (std::size_t index = 0; index < sizeof(names) / sizeof(names[0]); ++index)
  {
    if (index) out << '|';
    out << field(table, row, names[index]);
  }
  return out.str();
}

std::string eventKey(const Table& table, const std::vector<std::string>& row)
{
  const char* names[] = {
      "pair_generation", "pair_ordinal", "delivered_event_ordinal",
      "run_valid", "run_number", "event_valid", "event_number",
      "module_name"};
  std::ostringstream out;
  for (std::size_t index = 0; index < sizeof(names) / sizeof(names[0]); ++index)
  {
    if (index) out << '|';
    out << field(table, row, names[index]);
  }
  return out.str();
}

std::string eventScopeKey(const Table& table,
                          const std::vector<std::string>& row)
{
  return field(table, row, "pair_generation") + '|' +
         field(table, row, "pair_ordinal") + '|' +
         field(table, row, "delivered_event_ordinal") + '|' +
         field(table, row, "module_name");
}

std::string eventOccurrenceKey(const Table& table,
                               const std::vector<std::string>& row,
                               const char* prefix)
{
  const std::string base(prefix);
  return field(table, row, base + "_generation") + '|' +
         field(table, row, base + "_manager_ordinal") + '|' +
         field(table, row, base + "_ordinal");
}

struct Candidate
{
  std::uint64_t sequence;
  std::string identity;
  std::string event_identity;
  int photon_pt_slice;
  int centrality_slice;
  bool canonical;
  std::string tight_tag;
  bool preselection;
  bool admitted;
  std::uint16_t valid_mask;
  std::vector<std::string> variable_names;
  std::vector<std::string> value_bits;
  std::vector<std::string> triggers;
};

struct CandidateValueRow
{
  std::string variable_name;
  std::string value_bits;
  bool valid;
};

struct Inventory
{
  std::string directory;
  std::string object_name;
  std::string object_class;
  std::string title;
  std::string variable;
  std::string trigger;
  std::string tag;
  std::string view_suffix;
  int photon_pt_slice;
  int centrality_slice;
  int nbins;
  std::string xmin_bits;
  std::string xmax_bits;
  std::string x_axis_title;
  std::string y_axis_title;
  bool sumw2;
  bool required;
  std::string comparator;
};

std::string inventoryKey(const std::string& variable,
                         const std::string& trigger,
                         const std::string& tag,
                         int photon_pt_slice,
                         int centrality_slice)
{
  std::ostringstream out;
  out << variable << '|' << trigger << '|' << tag << '|'
      << photon_pt_slice << '|' << centrality_slice;
  return out.str();
}

std::string witnessKey(const std::string& candidate,
                       const std::string& variable,
                       const std::string& trigger,
                       const std::string& tag,
                       const std::string& object,
                       const std::string& value_bits,
                       const std::string& weight_bits,
                       int photon_pt_slice,
                       int centrality_slice,
                       bool sumw2)
{
  return candidate + '|' + variable + '|' + trigger + '|' + tag + '|' +
         object + '|' + value_bits + '|' + weight_bits + '|' +
         std::to_string(photon_pt_slice) + '|' +
         std::to_string(centrality_slice) + '|' + (sumw2 ? "1" : "0");
}

std::vector<Candidate> readCandidates(const std::string& cache_directory)
{
  const Table candidates = readTable(cache_directory + "/rawqa_candidates.tsv");
  const Table values = readTable(cache_directory + "/rawqa_candidate_values.tsv");
  const Table triggers = readTable(cache_directory + "/rawqa_candidate_triggers.tsv");
  std::map<std::uint64_t, std::map<std::size_t, CandidateValueRow> > value_rows;
  std::map<std::uint64_t, std::map<std::size_t, std::string> > trigger_rows;
  for (std::vector<std::vector<std::string> >::const_iterator row = values.rows.begin();
       row != values.rows.end(); ++row)
  {
    const std::uint64_t sequence = unsignedValue(field(values, *row, "candidate_sequence"));
    const std::size_t index = static_cast<std::size_t>(
        unsignedValue(field(values, *row, "value_index")));
    CandidateValueRow payload;
    payload.variable_name =
        the106::c0r::percentDecode(field(values, *row, "variable_name"));
    payload.value_bits = field(values, *row, "value_bits");
    payload.valid = booleanValue(field(values, *row, "valid"));
    if (!value_rows[sequence].insert(std::make_pair(index, payload)).second)
      throw std::runtime_error("duplicate candidate value index");
  }
  for (std::vector<std::vector<std::string> >::const_iterator row = triggers.rows.begin();
       row != triggers.rows.end(); ++row)
  {
    const std::uint64_t sequence = unsignedValue(field(triggers, *row, "candidate_sequence"));
    const std::size_t index = static_cast<std::size_t>(
        unsignedValue(field(triggers, *row, "trigger_index")));
    const std::string trigger =
        the106::c0r::percentDecode(field(triggers, *row, "trigger_name"));
    if (!trigger_rows[sequence].insert(std::make_pair(index, trigger)).second)
      throw std::runtime_error("duplicate candidate trigger index");
  }

  std::vector<Candidate> result;
  std::set<std::string> identities;
  std::set<std::uint64_t> candidate_sequences;
  for (std::vector<std::vector<std::string> >::const_iterator row = candidates.rows.begin();
       row != candidates.rows.end(); ++row)
  {
    Candidate candidate;
    candidate.sequence = unsignedValue(field(candidates, *row, "sequence"));
    candidate.identity = contextKey(candidates, *row);
    candidate.event_identity = eventKey(candidates, *row);
    candidate.photon_pt_slice = integerValue(field(candidates, *row, "photon_pt_slice"));
    candidate.centrality_slice = integerValue(field(candidates, *row, "centrality_slice"));
    candidate.canonical = booleanValue(field(candidates, *row, "canonical_view"));
    candidate.tight_tag = field(candidates, *row, "tight_tag");
    candidate.preselection = booleanValue(field(candidates, *row, "preselection_pass"));
    candidate.admitted = booleanValue(field(candidates, *row, "direct_candidate_admitted"));
    candidate.valid_mask = static_cast<std::uint16_t>(
        unsignedValue(field(candidates, *row, "valid_mask")));
    const std::size_t value_count = static_cast<std::size_t>(
        unsignedValue(field(candidates, *row, "value_count")));
    const std::size_t trigger_count = static_cast<std::size_t>(
        unsignedValue(field(candidates, *row, "active_trigger_count")));
    if (!identities.insert(candidate.identity).second)
      throw std::runtime_error("duplicate candidate identity");
    if (!candidate_sequences.insert(candidate.sequence).second)
      throw std::runtime_error("duplicate candidate sequence");
    for (std::size_t index = 0; index < value_count; ++index)
    {
      const std::map<std::size_t, CandidateValueRow>::const_iterator value =
          value_rows[candidate.sequence].find(index);
      if (value == value_rows[candidate.sequence].end())
        throw std::runtime_error("missing candidate value index");
      const bool mask_valid =
          (candidate.valid_mask & static_cast<std::uint16_t>(1U << index)) != 0;
      if (value->second.valid != mask_valid)
        throw std::runtime_error("candidate validity row disagrees with valid mask");
      candidate.variable_names.push_back(value->second.variable_name);
      candidate.value_bits.push_back(value->second.value_bits);
    }
    if (value_rows[candidate.sequence].size() != value_count)
      throw std::runtime_error("unexpected candidate value index");
    static const char* expected_variables[] = {
        "weta", "wphi", "weta33", "wphi33", "weta35",
        "wphi53", "et1", "e11e33", "e32e35"};
    if (value_count != sizeof(expected_variables) / sizeof(expected_variables[0]))
      throw std::runtime_error("candidate does not contain the exact nine-field domain");
    for (std::size_t index = 0; index < value_count; ++index)
      if (candidate.variable_names[index] != expected_variables[index])
        throw std::runtime_error("candidate raw-QA field order mismatch");
    for (std::size_t index = 0; index < trigger_count; ++index)
    {
      const std::map<std::size_t, std::string>::const_iterator trigger =
          trigger_rows[candidate.sequence].find(index);
      if (trigger == trigger_rows[candidate.sequence].end())
        throw std::runtime_error("missing candidate trigger index");
      candidate.triggers.push_back(trigger->second);
    }
    if (std::set<std::string>(candidate.triggers.begin(), candidate.triggers.end()).size() !=
        candidate.triggers.size())
      throw std::runtime_error("duplicate active trigger in candidate payload");
    if (trigger_rows[candidate.sequence].size() != trigger_count)
      throw std::runtime_error("unexpected candidate trigger index");
    result.push_back(candidate);
  }
  if (value_rows.size() != result.size())
    throw std::runtime_error("orphan candidate value rows");
  for (std::map<std::uint64_t, std::map<std::size_t, std::string> >::const_iterator
           row = trigger_rows.begin(); row != trigger_rows.end(); ++row)
    if (candidate_sequences.find(row->first) == candidate_sequences.end())
      throw std::runtime_error("orphan candidate trigger rows");
  return result;
}

std::map<std::string, Inventory> readInventory(const std::string& path)
{
  const Table table = readTable(path);
  std::map<std::string, Inventory> inventory;
  std::set<std::string> root_objects;
  for (std::vector<std::vector<std::string> >::const_iterator row = table.rows.begin();
       row != table.rows.end(); ++row)
  {
    Inventory item;
    item.directory = the106::c0r::percentDecode(field(table, *row, "directory"));
    item.object_name = the106::c0r::percentDecode(field(table, *row, "object_name"));
    item.object_class = the106::c0r::percentDecode(field(table, *row, "object_class"));
    item.title = the106::c0r::percentDecode(field(table, *row, "title"));
    item.variable = the106::c0r::percentDecode(field(table, *row, "variable"));
    item.trigger = the106::c0r::percentDecode(field(table, *row, "trigger"));
    item.tag = the106::c0r::percentDecode(field(table, *row, "tag"));
    item.view_suffix = the106::c0r::percentDecode(field(table, *row, "view_suffix"));
    item.photon_pt_slice = integerValue(field(table, *row, "photon_pt_slice"));
    item.centrality_slice = integerValue(field(table, *row, "centrality_slice"));
    item.nbins = integerValue(field(table, *row, "nbins"));
    item.xmin_bits = field(table, *row, "xmin_bits");
    item.xmax_bits = field(table, *row, "xmax_bits");
    item.x_axis_title = the106::c0r::percentDecode(field(table, *row, "x_axis_title"));
    item.y_axis_title = the106::c0r::percentDecode(field(table, *row, "y_axis_title"));
    item.sumw2 = booleanValue(field(table, *row, "sumw2_required"));
    item.required = booleanValue(field(table, *row, "required"));
    item.comparator = field(table, *row, "content_comparator");
    if (item.object_class != "TH1F" || item.nbins <= 0 ||
        item.directory.empty() || item.object_name.empty() || item.title.empty() ||
        item.variable.empty() || item.trigger.empty() || item.tag.empty() ||
        item.view_suffix != "canonical" || !item.required ||
        item.comparator != "BITWISE_SINGLE_PROCESS_UNWEIGHTED")
      throw std::runtime_error("unsupported C0-R inventory row");
    const std::string key = inventoryKey(
        item.variable, item.trigger, item.tag,
        item.photon_pt_slice, item.centrality_slice);
    if (!inventory.insert(std::make_pair(key, item)).second)
      throw std::runtime_error("duplicate inventory semantic key");
    if (!root_objects.insert(item.directory + "/" + item.object_name).second)
      throw std::runtime_error("duplicate inventory ROOT object");
  }
  if (inventory.empty())
    throw std::runtime_error("empty C0-R expected-object inventory");
  return inventory;
}

std::multiset<std::string> readWitnesses(const std::string& cache_directory)
{
  const Table table = readTable(cache_directory + "/rawqa_fill_witnesses.tsv");
  std::multiset<std::string> witnesses;
  std::map<std::string, std::uint64_t> expected_ordinal_by_event;
  for (std::vector<std::vector<std::string> >::const_iterator row = table.rows.begin();
       row != table.rows.end(); ++row)
  {
    if (!booleanValue(field(table, *row, "filled")) ||
        !booleanValue(field(table, *row, "value_valid")) ||
        !booleanValue(field(table, *row, "canonical_view")))
      throw std::runtime_error("invalid accepted fill witness");
    const std::string event_identity = eventKey(table, *row);
    const std::uint64_t ordinal = unsignedValue(field(table, *row, "fill_ordinal"));
    const std::uint64_t expected = ++expected_ordinal_by_event[event_identity];
    if (ordinal != expected)
      throw std::runtime_error("noncontiguous event-local fill-witness ordinal");
    const std::string candidate = contextKey(table, *row);
    const std::string directory =
        the106::c0r::percentDecode(field(table, *row, "directory_name"));
    const std::string object =
        the106::c0r::percentDecode(field(table, *row, "object_name"));
    const std::string object_class =
        the106::c0r::percentDecode(field(table, *row, "object_class"));
    const std::string contract =
        the106::c0r::percentDecode(field(table, *row, "object_contract_key"));
    if (object_class != "TH1F" || contract != directory + "/" + object)
      throw std::runtime_error("fill witness ROOT object contract disagreement");
    const std::string trigger =
        the106::c0r::percentDecode(field(table, *row, "trigger_name"));
    const std::string tag =
        the106::c0r::percentDecode(field(table, *row, "tag_name"));
    const std::string tight_tag = field(table, *row, "tight_tag");
    if (trigger != directory ||
        the106::c0r::percentDecode(field(table, *row, "view_suffix")) != "canonical" ||
        !booleanValue(field(table, *row, "canonical_view_context")) ||
        the106::c0r::percentDecode(field(table, *row, "view_key")) != "canonical")
      throw std::runtime_error("fill witness trigger/view contract disagreement");
    if ((tag == "tight" && tight_tag != "RawQATightTag::tight") ||
        (tag == "nonTight" && tight_tag != "RawQATightTag::non_tight") ||
        (tag != "pre" && tag != "tight" && tag != "nonTight"))
      throw std::runtime_error("fill witness tag/tight-state disagreement");
    const std::string key = witnessKey(
        candidate,
        the106::c0r::percentDecode(field(table, *row, "variable_name")),
        trigger,
        tag,
        object,
        field(table, *row, "value_bits"),
        field(table, *row, "weight_bits"),
        integerValue(field(table, *row, "photon_pt_slice")),
        integerValue(field(table, *row, "centrality_slice")),
        booleanValue(field(table, *row, "sumw2_enabled")));
    witnesses.insert(key);
  }
  return witnesses;
}

void validateCompletion(const std::string& cache_directory)
{
  const Table metadata = readTable(cache_directory + "/metadata.tsv");
  const Table completion = readTable(cache_directory + "/completion.tsv");
  std::map<std::string, std::string> metadata_values;
  std::map<std::string, std::string> completion_values;
  for (std::vector<std::vector<std::string> >::const_iterator row = metadata.rows.begin();
       row != metadata.rows.end(); ++row)
    metadata_values[field(metadata, *row, "key")] =
        the106::c0r::percentDecode(field(metadata, *row, "value"));
  for (std::vector<std::vector<std::string> >::const_iterator row = completion.rows.begin();
       row != completion.rows.end(); ++row)
    completion_values[field(completion, *row, "key")] = field(completion, *row, "value");
  if (metadata_values["format"] != "THE106_C0R_CACHE_V1")
    throw std::runtime_error("wrong cache format");
  const char* required_metadata[] = {
      "canary_id", "scientific_contract_sha256", "source_manifest_sha256",
      "configuration_sha256", "reconstruction_identity", "cdb_identity",
      "model_reference", "model_sha256", "feature_order_sha256",
      "preprocessing_identity", "runtime_provider_identity", "source_range"};
  for (std::size_t index = 0;
       index < sizeof(required_metadata) / sizeof(required_metadata[0]); ++index)
    if (metadata_values[required_metadata[index]].empty())
      throw std::runtime_error(std::string("missing sealed metadata: ") +
                               required_metadata[index]);
  const std::string source_range = metadata_values["source_range"];
  if (source_range.size() < 5 || source_range.compare(0, 3, "[0,") != 0 ||
      source_range[source_range.size() - 1] != ']')
    throw std::runtime_error("cache source range is not canonical [0,N]");
  const char* required_zero[] = {
      "source_observer_failures", "source_invariant_violations",
      "recoil_observer_failures", "recoil_invariant_violations",
      "serialization_failures"};
  for (std::size_t index = 0; index < sizeof(required_zero) / sizeof(required_zero[0]); ++index)
  {
    if (completion_values[required_zero[index]] != "0")
      throw std::runtime_error(std::string("nonzero completion failure: ") + required_zero[index]);
  }
  if (completion_values["closed_world"] != "1")
    throw std::runtime_error("cache source population is not closed world");
}

std::uint64_t expectedSourcePairCount(const std::string& cache_directory)
{
  const Table metadata = readTable(cache_directory + "/metadata.tsv");
  std::string source_range;
  for (std::vector<std::vector<std::string> >::const_iterator row =
           metadata.rows.begin(); row != metadata.rows.end(); ++row)
    if (field(metadata, *row, "key") == "source_range")
      source_range = the106::c0r::percentDecode(field(metadata, *row, "value"));
  if (source_range.size() < 5 || source_range.compare(0, 3, "[0,") != 0 ||
      source_range[source_range.size() - 1] != ']')
    throw std::runtime_error("missing canonical source range");
  const std::uint64_t upper = unsignedValue(
      source_range.substr(3, source_range.size() - 4));
  if (upper == std::numeric_limits<std::uint64_t>::max())
    throw std::runtime_error("source range overflows pair count");
  return upper + 1;
}

std::string occurrenceKey(const Table& table, const std::vector<std::string>& row)
{
  return field(table, row, "occurrence_generation") + '|' +
         field(table, row, "manager_ordinal") + '|' +
         field(table, row, "occurrence_ordinal");
}

std::string pairKey(const Table& table, const std::vector<std::string>& row)
{
  return field(table, row, "pair_generation") + '|' +
         field(table, row, "pair_ordinal");
}

struct SourceClosure
{
  std::set<std::string> delivered_pairs;
  std::map<std::string, std::set<std::string> > delivered_occurrences;
};

SourceClosure validateSource(const std::string& cache_directory,
                             std::uint64_t expected_pair_count)
{
  const Table table = readTable(cache_directory + "/source_observations.tsv");
  std::set<std::string> occurrence_begins;
  std::set<std::string> occurrence_terminals;
  std::set<std::string> pair_begins;
  std::set<std::string> pair_terminals;
  std::map<std::string, std::set<std::string> > delivered_occurrences;
  std::map<std::uint64_t, std::set<std::uint64_t> > occurrence_ordinals;
  std::map<std::uint64_t, std::set<std::uint64_t> > source_entries;
  std::map<std::uint64_t, std::set<std::string> > source_descriptors;
  std::map<std::string, std::set<std::uint64_t> > read_attempts;
  std::map<std::string, std::string> pair_terminal_disposition;
  std::set<std::uint64_t> pair_ordinals;
  for (std::vector<std::vector<std::string> >::const_iterator row = table.rows.begin();
       row != table.rows.end(); ++row)
  {
    const std::string kind = field(table, *row, "kind");
    const std::string occurrence = occurrenceKey(table, *row);
    const std::string pair = pairKey(table, *row);
    if (kind == "SourceEventKind::occurrence_begin")
    {
      if (!occurrence_begins.insert(occurrence).second)
        throw std::runtime_error("duplicate logical occurrence begin");
      occurrence_ordinals[unsignedValue(field(table, *row, "manager_ordinal"))]
          .insert(unsignedValue(field(table, *row, "occurrence_ordinal")));
      source_entries[unsignedValue(field(table, *row, "manager_ordinal"))]
          .insert(unsignedValue(field(table, *row, "source_entry")));
      source_descriptors[unsignedValue(field(table, *row, "manager_ordinal"))]
          .insert(the106::c0r::percentDecode(
              field(table, *row, "source_descriptor")));
    }
    else if (kind == "SourceEventKind::occurrence_resume" ||
             kind == "SourceEventKind::low_level_read" ||
             kind == "SourceEventKind::occurrence_transition")
    {
      if (occurrence_begins.find(occurrence) == occurrence_begins.end())
        throw std::runtime_error("read/resume/transition without occurrence begin");
      if (kind == "SourceEventKind::low_level_read")
      {
        const std::uint64_t attempt = unsignedValue(
            field(table, *row, "read_attempt_ordinal"));
        if (attempt == 0 || !read_attempts[occurrence].insert(attempt).second)
          throw std::runtime_error("duplicate or zero low-level read attempt");
      }
    }
    else if (kind == "SourceEventKind::occurrence_terminal")
    {
      if (occurrence_begins.find(occurrence) == occurrence_begins.end() ||
          !occurrence_terminals.insert(occurrence).second)
        throw std::runtime_error("missing or duplicate occurrence terminal");
      if (field(table, *row, "terminal") ==
          "OccurrenceTerminal::delivered_to_modules")
      {
        if (field(table, *row, "pair_generation") == "0" ||
            field(table, *row, "pair_ordinal") == "0")
          throw std::runtime_error("delivered occurrence lacks pair token");
        delivered_occurrences[pair].insert(occurrence);
      }
    }
    else if (kind == "SourceEventKind::pair_begin")
    {
      if (!pair_begins.insert(pair).second)
        throw std::runtime_error("duplicate pair begin");
      pair_ordinals.insert(unsignedValue(field(table, *row, "pair_ordinal")));
    }
    else if (kind == "SourceEventKind::pair_transition")
    {
      if (pair_begins.find(pair) == pair_begins.end())
        throw std::runtime_error("pair transition without pair begin");
    }
    else if (kind == "SourceEventKind::pair_terminal")
    {
      if (pair_begins.find(pair) == pair_begins.end() ||
          !pair_terminals.insert(pair).second)
        throw std::runtime_error("missing or duplicate pair terminal");
      pair_terminal_disposition[pair] = field(table, *row, "terminal");
    }
    else if (kind == "SourceEventKind::invariant_violation")
      throw std::runtime_error("source observer emitted invariant violation");
    else
      throw std::runtime_error("unknown source observation kind");
  }
  if (occurrence_begins != occurrence_terminals)
    throw std::runtime_error("source occurrence population is not closed world");
  if (pair_begins != pair_terminals)
    throw std::runtime_error("source pair population is not closed world");
  if (occurrence_ordinals.size() != 2)
    throw std::runtime_error("frozen paired source does not contain exactly two managers");
  std::set<std::uint64_t> all_occurrence_ordinals;
  for (std::map<std::uint64_t, std::set<std::uint64_t> >::const_iterator manager =
           occurrence_ordinals.begin(); manager != occurrence_ordinals.end(); ++manager)
  {
    if (manager->first == 0 || manager->second.size() != expected_pair_count)
      throw std::runtime_error("manager occurrence count disagrees with frozen range");
    all_occurrence_ordinals.insert(manager->second.begin(), manager->second.end());
    if (source_entries[manager->first].size() != expected_pair_count)
      throw std::runtime_error("source-entry count disagrees with frozen range");
    for (std::uint64_t entry = 0; entry < expected_pair_count; ++entry)
      if (source_entries[manager->first].find(entry) ==
          source_entries[manager->first].end())
        throw std::runtime_error("source-entry range is not contiguous [0,N]");
    if (source_descriptors[manager->first].size() != 1 ||
        source_descriptors[manager->first].begin()->empty())
      throw std::runtime_error("manager lacks one immutable source descriptor");
  }
  if (all_occurrence_ordinals.size() != 2 * expected_pair_count)
    throw std::runtime_error("logical occurrence ordinals collide across source managers");
  for (std::uint64_t ordinal = 1; ordinal <= 2 * expected_pair_count; ++ordinal)
    if (all_occurrence_ordinals.find(ordinal) == all_occurrence_ordinals.end())
      throw std::runtime_error("global logical occurrence range is not contiguous");
  if (pair_ordinals.size() != expected_pair_count)
    throw std::runtime_error("pair count disagrees with frozen range");
  for (std::uint64_t ordinal = 1; ordinal <= expected_pair_count; ++ordinal)
    if (pair_ordinals.find(ordinal) == pair_ordinals.end())
      throw std::runtime_error("pair ordinal range is not contiguous");
  for (std::set<std::string>::const_iterator occurrence = occurrence_begins.begin();
       occurrence != occurrence_begins.end(); ++occurrence)
  {
    const std::set<std::uint64_t>& attempts = read_attempts[*occurrence];
    if (attempts.empty())
      throw std::runtime_error("logical occurrence has no low-level read observation");
    const std::uint64_t last = *attempts.rbegin();
    for (std::uint64_t attempt = 1; attempt <= last; ++attempt)
      if (attempts.find(attempt) == attempts.end())
        throw std::runtime_error("low-level read-attempt range is not contiguous");
  }
  SourceClosure closure;
  for (std::map<std::string, std::set<std::string> >::const_iterator pair =
           delivered_occurrences.begin(); pair != delivered_occurrences.end(); ++pair)
  {
    if (pair->second.size() != 2)
      throw std::runtime_error("delivered pair does not contain exactly two stream occurrences");
    if (pair_terminal_disposition[pair->first] !=
        "OccurrenceTerminal::delivered_to_modules")
      throw std::runtime_error("delivered occurrence pair has non-delivery terminal");
    closure.delivered_pairs.insert(pair->first);
    closure.delivered_occurrences.insert(*pair);
  }
  return closure;
}

void validateEvents(const std::string& cache_directory,
                    const std::vector<Candidate>& candidates,
                    const SourceClosure& source)
{
  const Table table = readTable(cache_directory + "/event_observations.tsv");
  const Table trigger_table = readTable(
      cache_directory + "/event_active_triggers.tsv");
  const Table fill_table = readTable(
      cache_directory + "/rawqa_fill_witnesses.tsv");
  std::map<std::uint64_t, std::map<std::size_t, std::string> > trigger_rows;
  for (std::vector<std::vector<std::string> >::const_iterator row =
           trigger_table.rows.begin(); row != trigger_table.rows.end(); ++row)
  {
    const std::uint64_t sequence = unsignedValue(
        field(trigger_table, *row, "event_sequence"));
    const std::size_t index = static_cast<std::size_t>(unsignedValue(
        field(trigger_table, *row, "trigger_index")));
    const std::string trigger = the106::c0r::percentDecode(
        field(trigger_table, *row, "trigger_name"));
    if (!trigger_rows[sequence].insert(std::make_pair(index, trigger)).second)
      throw std::runtime_error("duplicate event active-trigger index");
  }
  std::map<std::string, std::vector<std::string> > finals;
  std::map<std::string, std::vector<std::string> > event_triggers;
  std::set<std::string> final_pairs;
  std::map<std::string, std::uint64_t> final_pair_counts;
  std::map<std::string, std::map<std::string, std::uint64_t> > scope_kinds;
  std::set<std::uint64_t> final_sequences;
  for (std::vector<std::vector<std::string> >::const_iterator row = table.rows.begin();
       row != table.rows.end(); ++row)
  {
    const std::string kind = field(table, *row, "kind");
    if (kind != "EventObservationKind::begin" &&
        kind != "EventObservationKind::context" &&
        kind != "EventObservationKind::finalization")
      throw std::runtime_error("unknown event-observation kind");
    ++scope_kinds[eventScopeKey(table, *row)][kind];
    if (kind != "EventObservationKind::finalization") continue;
    const std::uint64_t sequence = unsignedValue(field(table, *row, "sequence"));
    final_sequences.insert(sequence);
    const std::string key = eventKey(table, *row);
    const std::string pair = field(table, *row, "pair_generation") + '|' +
                             field(table, *row, "pair_ordinal");
    if (source.delivered_pairs.find(pair) == source.delivered_pairs.end())
      throw std::runtime_error("event finalization lacks a delivered source pair");
    final_pairs.insert(pair);
    ++final_pair_counts[pair];
    if (unsignedValue(field(table, *row, "occurrence_count")) != 2)
      throw std::runtime_error("event finalization does not bind two source occurrences");
    std::set<std::string> event_occurrences;
    event_occurrences.insert(eventOccurrenceKey(table, *row, "occurrence0"));
    event_occurrences.insert(eventOccurrenceKey(table, *row, "occurrence1"));
    const std::map<std::string, std::set<std::string> >::const_iterator expected_occurrences =
        source.delivered_occurrences.find(pair);
    if (expected_occurrences == source.delivered_occurrences.end() ||
        event_occurrences != expected_occurrences->second)
      throw std::runtime_error("event finalization source-occurrence linkage mismatch");
    if (!finals.insert(std::make_pair(key, *row)).second)
      throw std::runtime_error("duplicate event finalization");
    if (!booleanValue(field(table, *row, "finalized")))
      throw std::runtime_error("non-final event finalization row");
    const std::uint64_t admitted = unsignedValue(field(table, *row, "admitted_candidate_count"));
    const std::uint64_t fills = unsignedValue(field(table, *row, "rawqa_fill_count"));
    const std::size_t trigger_count = static_cast<std::size_t>(unsignedValue(
        field(table, *row, "active_trigger_count")));
    for (std::size_t index = 0; index < trigger_count; ++index)
    {
      if (trigger_rows[sequence].find(index) == trigger_rows[sequence].end())
        throw std::runtime_error("event finalization lacks active-trigger row");
      event_triggers[key].push_back(trigger_rows[sequence].find(index)->second);
    }
    if (trigger_rows[sequence].size() != trigger_count)
      throw std::runtime_error("event finalization has unexpected active-trigger row");
    if (booleanValue(field(table, *row, "zero_candidate")) != (admitted == 0) ||
        booleanValue(field(table, *row, "zero_fill")) != (fills == 0))
      throw std::runtime_error("explicit zero-candidate/fill state disagreement");
  }
  if (finals.empty()) throw std::runtime_error("no event finalizations");
  if (final_pairs != source.delivered_pairs)
    throw std::runtime_error("delivered source pairs and event finalizations differ");
  for (std::map<std::string, std::uint64_t>::const_iterator pair =
           final_pair_counts.begin(); pair != final_pair_counts.end(); ++pair)
    if (pair->second != 1)
      throw std::runtime_error("delivered pair has other than one event finalization");
  for (std::map<std::string, std::map<std::string, std::uint64_t> >::const_iterator scope =
           scope_kinds.begin(); scope != scope_kinds.end(); ++scope)
  {
    if (scope->second.find("EventObservationKind::begin") == scope->second.end() ||
        scope->second.find("EventObservationKind::begin")->second != 1 ||
        scope->second.find("EventObservationKind::finalization") == scope->second.end() ||
        scope->second.find("EventObservationKind::finalization")->second != 1)
      throw std::runtime_error("event scope lacks exactly one begin and finalization");
    const std::map<std::string, std::uint64_t>::const_iterator contexts =
        scope->second.find("EventObservationKind::context");
    if (contexts != scope->second.end() && contexts->second > 1)
      throw std::runtime_error("event scope has duplicate context observations");
  }
  for (std::map<std::uint64_t, std::map<std::size_t, std::string> >::const_iterator
           row = trigger_rows.begin(); row != trigger_rows.end(); ++row)
    if (final_sequences.find(row->first) == final_sequences.end())
      throw std::runtime_error("active-trigger rows do not belong to a finalization");
  std::map<std::string, std::uint64_t> observed_candidates;
  std::set<std::string> candidate_identities;
  for (std::vector<Candidate>::const_iterator candidate = candidates.begin();
       candidate != candidates.end(); ++candidate)
  {
    if (finals.find(candidate->event_identity) == finals.end())
      throw std::runtime_error("candidate lacks exact event finalization");
    ++observed_candidates[candidate->event_identity];
    candidate_identities.insert(candidate->identity);
    if (candidate->triggers != event_triggers[candidate->event_identity])
      throw std::runtime_error("candidate and event active-trigger payloads differ");
  }
  std::map<std::string, std::uint64_t> observed_fills;
  std::map<std::string, std::set<std::string> > observed_contributors;
  for (std::vector<std::vector<std::string> >::const_iterator row =
           fill_table.rows.begin(); row != fill_table.rows.end(); ++row)
  {
    if (!booleanValue(field(fill_table, *row, "filled"))) continue;
    const std::string event_identity = eventKey(fill_table, *row);
    const std::string candidate_identity = contextKey(fill_table, *row);
    if (finals.find(event_identity) == finals.end())
      throw std::runtime_error("fill witness lacks exact event finalization");
    if (candidate_identities.find(candidate_identity) == candidate_identities.end())
      throw std::runtime_error("fill witness lacks exact candidate row");
    ++observed_fills[event_identity];
    observed_contributors[event_identity].insert(candidate_identity);
  }
  for (std::map<std::string, std::vector<std::string> >::const_iterator
           event = finals.begin(); event != finals.end(); ++event)
  {
    const std::uint64_t admitted = unsignedValue(
        field(table, event->second, "admitted_candidate_count"));
    if (admitted != observed_candidates[event->first])
      throw std::runtime_error("event admitted-candidate count disagrees with cache rows");
    const std::uint64_t fill_count = unsignedValue(
        field(table, event->second, "rawqa_fill_count"));
    const std::uint64_t contributor_count = unsignedValue(
        field(table, event->second, "rawqa_contributor_count"));
    if (fill_count != observed_fills[event->first] ||
        contributor_count != observed_contributors[event->first].size())
      throw std::runtime_error("event fill/contributor counts disagree with witnesses");
  }
}

struct ReplaySummary
{
  std::uint64_t candidate_count;
  std::uint64_t predicted_fill_count;
  std::uint64_t witness_count;
  std::uint64_t object_count;
};

ReplaySummary replay(const std::string& cache_directory,
                     const std::string& inventory_path,
                     const std::string& output_path)
{
  validateCompletion(cache_directory);
  const SourceClosure source = validateSource(
      cache_directory, expectedSourcePairCount(cache_directory));
  const std::vector<Candidate> candidates = readCandidates(cache_directory);
  validateEvents(cache_directory, candidates, source);
  const std::map<std::string, Inventory> inventory = readInventory(inventory_path);
  std::multiset<std::string> witnesses = readWitnesses(cache_directory);

  std::unique_ptr<TFile> output(TFile::Open(output_path.c_str(), "CREATE"));
  if (!output || output->IsZombie())
    throw std::runtime_error("refusing to overwrite or unable to create replay ROOT file");

  std::map<std::string, TH1F*> histograms;
  for (std::map<std::string, Inventory>::const_iterator iter = inventory.begin();
       iter != inventory.end(); ++iter)
  {
    const Inventory& item = iter->second;
    double xmin = 0.0;
    double xmax = 0.0;
    if (!the106::c0r::binary64FromHex(item.xmin_bits, xmin) ||
        !the106::c0r::binary64FromHex(item.xmax_bits, xmax) || !(xmin < xmax))
      throw std::runtime_error("invalid inventory binary64 axis boundary");
    TDirectory* directory = output->GetDirectory(item.directory.c_str());
    if (!directory) directory = output->mkdir(item.directory.c_str());
    if (!directory) throw std::runtime_error("unable to create ROOT directory");
    directory->cd();
    TH1F* histogram = new TH1F(
        item.object_name.c_str(), item.title.c_str(), item.nbins, xmin, xmax);
    histogram->GetXaxis()->SetTitle(item.x_axis_title.c_str());
    histogram->GetYaxis()->SetTitle(item.y_axis_title.c_str());
    if (item.sumw2) histogram->Sumw2();
    histogram->SetDirectory(directory);
    histograms[iter->first] = histogram;
    output->cd();
  }

  const std::string unit_weight_bits = the106::c0r::binary64Hex(1.0);
  std::uint64_t predicted_fills = 0;
  for (std::vector<Candidate>::const_iterator candidate = candidates.begin();
       candidate != candidates.end(); ++candidate)
  {
    if (!candidate->canonical || !candidate->preselection || !candidate->admitted)
      throw std::runtime_error("candidate outside the frozen canonical admitted domain");
    std::vector<std::string> tags(1, "pre");
    if (candidate->tight_tag == "RawQATightTag::tight") tags.push_back("tight");
    else if (candidate->tight_tag == "RawQATightTag::non_tight") tags.push_back("nonTight");
    else if (candidate->tight_tag != "RawQATightTag::neither")
      throw std::runtime_error("unsupported raw-QA tight tag");

    if (candidate->variable_names.size() != 9 || candidate->value_bits.size() != 9)
      throw std::runtime_error("raw-QA candidate does not have nine values");
    for (std::size_t trigger_index = 0;
         trigger_index < candidate->triggers.size(); ++trigger_index)
    {
      for (std::size_t tag_index = 0; tag_index < tags.size(); ++tag_index)
      {
        for (std::size_t value_index = 0; value_index < 9; ++value_index)
        {
          if ((candidate->valid_mask & (1U << value_index)) == 0) continue;
          double value = 0.0;
          if (!the106::c0r::binary64FromHex(candidate->value_bits[value_index], value) ||
              !std::isfinite(value))
            throw std::runtime_error("valid raw-QA value has invalid binary64 payload");
          const std::string semantic_key = inventoryKey(
              candidate->variable_names[value_index], candidate->triggers[trigger_index],
              tags[tag_index], candidate->photon_pt_slice,
              candidate->centrality_slice);
          const std::map<std::string, Inventory>::const_iterator contract =
              inventory.find(semantic_key);
          if (contract == inventory.end())
            throw std::runtime_error("candidate contribution absent from frozen inventory");
          const std::string witness = witnessKey(
              candidate->identity, candidate->variable_names[value_index],
              candidate->triggers[trigger_index], tags[tag_index],
              contract->second.object_name, candidate->value_bits[value_index],
              unit_weight_bits, candidate->photon_pt_slice,
              candidate->centrality_slice, contract->second.sumw2);
          const std::multiset<std::string>::iterator found = witnesses.find(witness);
          if (found == witnesses.end())
            throw std::runtime_error("candidate-derived contribution lacks exact fill witness");
          witnesses.erase(found);
          histograms.at(semantic_key)->Fill(value);
          ++predicted_fills;
        }
      }
    }
  }
  if (!witnesses.empty())
    throw std::runtime_error("fill witnesses contain non-candidate-derived contributions");

  output->Write();
  output->Close();
  ReplaySummary summary;
  summary.candidate_count = candidates.size();
  summary.predicted_fill_count = predicted_fills;
  summary.witness_count = predicted_fills;
  summary.object_count = inventory.size();
  return summary;
}

struct CompareSummary
{
  std::uint64_t object_count;
  std::uint64_t bin_count;
};

CompareSummary compare(const std::string& inventory_path,
                       const std::string& lhs_path,
                       const std::string& rhs_path)
{
  const std::map<std::string, Inventory> inventory = readInventory(inventory_path);
  std::unique_ptr<TFile> lhs(TFile::Open(lhs_path.c_str(), "READ"));
  std::unique_ptr<TFile> rhs(TFile::Open(rhs_path.c_str(), "READ"));
  if (!lhs || lhs->IsZombie() || !rhs || rhs->IsZombie())
    throw std::runtime_error("unable to open ROOT files for closure comparison");
  std::map<std::string, std::set<std::string> > expected_names_by_directory;
  std::set<std::string> expected_root_paths;
  for (std::map<std::string, Inventory>::const_iterator iter = inventory.begin();
       iter != inventory.end(); ++iter)
  {
    expected_names_by_directory[iter->second.directory].insert(iter->second.object_name);
    expected_root_paths.insert(iter->second.directory + "/" + iter->second.object_name);
  }

  const char* variables[] = {
      "weta", "wphi", "weta33", "wphi33", "weta35",
      "wphi53", "et1", "e11e33", "e32e35"};
  const char* tags[] = {"pre", "tight", "nonTight"};
  const auto is_target_name = [&](const std::string& name) -> bool
  {
    for (std::size_t variable = 0;
         variable < sizeof(variables) / sizeof(variables[0]); ++variable)
      for (std::size_t tag = 0; tag < sizeof(tags) / sizeof(tags[0]); ++tag)
        if (name.find(std::string("h_ss_") + variables[variable] + "_" +
                      tags[tag] + "_") == 0)
          return true;
    return false;
  };
  TFile* files[] = {lhs.get(), rhs.get()};
  for (std::size_t side = 0; side < 2; ++side)
  {
    TIter next_top(files[side]->GetListOfKeys());
    while (TKey* top_key = dynamic_cast<TKey*>(next_top()))
    {
      TDirectory* directory = files[side]->GetDirectory(top_key->GetName());
      if (!directory) continue;
      TIter next_object(directory->GetListOfKeys());
      while (TKey* object_key = dynamic_cast<TKey*>(next_object()))
      {
        const std::string name = object_key->GetName();
        if (!is_target_name(name)) continue;
        const std::string path = std::string(top_key->GetName()) + "/" + name;
        if (expected_root_paths.find(path) == expected_root_paths.end())
          throw std::runtime_error("frozen inventory omits target ROOT object " + path);
      }
    }
  }
  for (std::map<std::string, std::set<std::string> >::const_iterator directory =
           expected_names_by_directory.begin();
       directory != expected_names_by_directory.end(); ++directory)
  {
    TDirectory* left_directory = lhs->GetDirectory(directory->first.c_str());
    TDirectory* right_directory = rhs->GetDirectory(directory->first.c_str());
    if (!left_directory || !right_directory)
      throw std::runtime_error("missing frozen ROOT directory " + directory->first);
    TDirectory* directories[] = {left_directory, right_directory};
    for (std::size_t side = 0; side < 2; ++side)
    {
      TIter next(directories[side]->GetListOfKeys());
      while (TKey* key = dynamic_cast<TKey*>(next()))
      {
        const std::string name = key->GetName();
        if (is_target_name(name) && directory->second.find(name) == directory->second.end())
          throw std::runtime_error("unexpected frozen-family ROOT object " +
                                   directory->first + "/" + name);
      }
    }
  }
  std::uint64_t bins = 0;
  for (std::map<std::string, Inventory>::const_iterator iter = inventory.begin();
       iter != inventory.end(); ++iter)
  {
    const Inventory& item = iter->second;
    const std::string path = item.directory + "/" + item.object_name;
    TH1F* left = dynamic_cast<TH1F*>(lhs->Get(path.c_str()));
    TH1F* right = dynamic_cast<TH1F*>(rhs->Get(path.c_str()));
    if (!left || !right) throw std::runtime_error("missing required ROOT object " + path);
    TDirectory* left_directory = lhs->GetDirectory(item.directory.c_str());
    TDirectory* right_directory = rhs->GetDirectory(item.directory.c_str());
    std::size_t left_key_count = 0;
    std::size_t right_key_count = 0;
    TIter next_left(left_directory->GetListOfKeys());
    while (TKey* key = dynamic_cast<TKey*>(next_left()))
      if (item.object_name == key->GetName()) ++left_key_count;
    TIter next_right(right_directory->GetListOfKeys());
    while (TKey* key = dynamic_cast<TKey*>(next_right()))
      if (item.object_name == key->GetName()) ++right_key_count;
    if (left_key_count != 1 || right_key_count != 1)
      throw std::runtime_error("required ROOT object does not appear exactly once " + path);
    double xmin = 0.0;
    double xmax = 0.0;
    if (!the106::c0r::binary64FromHex(item.xmin_bits, xmin) ||
        !the106::c0r::binary64FromHex(item.xmax_bits, xmax))
      throw std::runtime_error("invalid frozen axis bits " + path);
    if (std::string(left->ClassName()) != item.object_class ||
        std::string(right->ClassName()) != item.object_class ||
        std::string(left->GetTitle()) != item.title ||
        std::string(right->GetTitle()) != item.title ||
        left->GetNbinsX() != item.nbins || right->GetNbinsX() != item.nbins ||
        left->GetSumw2N() != right->GetSumw2N() ||
        left->GetEntries() != right->GetEntries() ||
        the106::c0r::binary64Hex(left->GetXaxis()->GetBinLowEdge(1)) !=
            the106::c0r::binary64Hex(xmin) ||
        the106::c0r::binary64Hex(right->GetXaxis()->GetBinLowEdge(1)) !=
            the106::c0r::binary64Hex(xmin) ||
        the106::c0r::binary64Hex(left->GetXaxis()->GetBinUpEdge(item.nbins)) !=
            the106::c0r::binary64Hex(xmax) ||
        the106::c0r::binary64Hex(right->GetXaxis()->GetBinUpEdge(item.nbins)) !=
            the106::c0r::binary64Hex(xmax) ||
        std::string(left->GetXaxis()->GetTitle()) != item.x_axis_title ||
        std::string(right->GetXaxis()->GetTitle()) != item.x_axis_title ||
        std::string(left->GetYaxis()->GetTitle()) != item.y_axis_title ||
        std::string(right->GetYaxis()->GetTitle()) != item.y_axis_title)
      throw std::runtime_error("ROOT object contract mismatch " + path);
    if (left->GetSumw2N() != (item.sumw2 ? item.nbins + 2 : 0) ||
        right->GetSumw2N() != (item.sumw2 ? item.nbins + 2 : 0))
      throw std::runtime_error("ROOT Sumw2 contract mismatch " + path);
    for (int bin = 0; bin <= item.nbins + 1; ++bin)
    {
      if (left->GetBinContent(bin) != right->GetBinContent(bin) ||
          left->GetBinError(bin) != right->GetBinError(bin))
        throw std::runtime_error("ROOT content/Sumw2 mismatch " + path);
      if (item.sumw2 &&
          left->GetSumw2()->At(bin) != right->GetSumw2()->At(bin))
        throw std::runtime_error("ROOT raw Sumw2 mismatch " + path);
      ++bins;
    }
  }
  CompareSummary summary;
  summary.object_count = inventory.size();
  summary.bin_count = bins;
  return summary;
}

void writeReplayReport(const std::string& path,
                       const ReplaySummary& summary)
{
  std::ofstream out(path.c_str());
  if (!out) throw std::runtime_error("unable to create replay report");
  out << "{\n"
      << "  \"result\": \"PASS\",\n"
      << "  \"route\": \"CACHE_ONLY_DIAGNOSTIC_REPLAY\",\n"
      << "  \"candidate_count\": " << summary.candidate_count << ",\n"
      << "  \"predicted_fill_count\": " << summary.predicted_fill_count << ",\n"
      << "  \"witness_count\": " << summary.witness_count << ",\n"
      << "  \"object_count\": " << summary.object_count << "\n"
      << "}\n";
}

void writeCompareReport(const std::string& path,
                        const CompareSummary& summary,
                        const std::string& lhs,
                        const std::string& rhs)
{
  std::ofstream out(path.c_str());
  if (!out) throw std::runtime_error("unable to create comparison report");
  out << "{\n"
      << "  \"result\": \"PASS\",\n"
      << "  \"comparison\": \"BITWISE_SINGLE_PROCESS_ROOT_CONTENT\",\n"
      << "  \"lhs\": \"" << jsonEscape(lhs) << "\",\n"
      << "  \"rhs\": \"" << jsonEscape(rhs) << "\",\n"
      << "  \"object_count\": " << summary.object_count << ",\n"
      << "  \"bin_count_including_flow\": " << summary.bin_count << "\n"
      << "}\n";
}

void usage(const char* program)
{
  std::cerr
      << "usage:\n"
      << "  " << program << " replay CACHE_DIR INVENTORY.tsv OUTPUT.root REPORT.json\n"
      << "  " << program << " compare INVENTORY.tsv LHS.root RHS.root REPORT.json\n";
}

}  // namespace

int main(int argc, char** argv)
{
  try
  {
    if (argc == 6 && std::string(argv[1]) == "replay")
    {
      const ReplaySummary summary = replay(argv[2], argv[3], argv[4]);
      writeReplayReport(argv[5], summary);
      return 0;
    }
    if (argc == 6 && std::string(argv[1]) == "compare")
    {
      const CompareSummary summary = compare(argv[2], argv[3], argv[4]);
      writeCompareReport(argv[5], summary, argv[3], argv[4]);
      return 0;
    }
    usage(argv[0]);
    return 2;
  }
  catch (const std::exception& error)
  {
    std::cerr << "THE106 C0-R fail-closed: " << error.what() << '\n';
    return 1;
  }
}
