#include "the106_c0r_diagnostic_writer.h"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <cstdlib>
#include <unistd.h>

namespace
{

std::size_t lineCount(const std::string& path)
{
  std::ifstream input(path.c_str());
  assert(input.good());
  std::size_t count = 0;
  std::string line;
  while (std::getline(input, line)) ++count;
  return count;
}

std::uint64_t bits(double value)
{
  std::uint64_t result = 0;
  std::memcpy(&result, &value, sizeof(result));
  return result;
}

}  // namespace

int main()
{
  char temporary[] = "/tmp/the106_c0r_writer_test_XXXXXX";
  const char* requested_directory = std::getenv("THE106_C0R_TEST_CACHE_DIR");
  const char* directory_value = requested_directory;
  if (!directory_value)
  {
    directory_value = ::mkdtemp(temporary);
    assert(directory_value);
  }
  const std::string directory(directory_value);

  the106::c0r::RunMetadata metadata;
  metadata.canary_id = "synthetic-C0-R";
  metadata.scientific_contract_sha256 = std::string(64, 'a');
  metadata.source_manifest_sha256 = std::string(64, 'b');
  metadata.configuration_sha256 = std::string(64, 'c');
  metadata.reconstruction_identity = "synthetic";
  metadata.cdb_identity = "synthetic-cdb";
  metadata.model_reference = "model.root";
  metadata.model_sha256 = std::string(64, 'd');
  metadata.feature_order_sha256 = std::string(64, 'e');
  metadata.preprocessing_identity = "identity";
  metadata.runtime_provider_identity = "TMVA::Experimental::RBDT";
  metadata.source_range = "[0,0]";

  {
    the106::c0r::DiagnosticWriter writer(directory, metadata);
    assert(writer.good());

    const void* sync = reinterpret_cast<const void*>(0x100);
    const void* manager_a = reinterpret_cast<const void*>(0x101);
    const void* manager_b = reinterpret_cast<const void*>(0x102);
    const void* io_a = reinterpret_cast<const void*>(0x201);
    const void* io_b = reinterpret_cast<const void*>(0x202);
    the106::c0h2::bindInputManagerIO(
        manager_a, io_a, sync, "DST_JET", "physical-source-a");
    the106::c0h2::bindInputManagerIO(
        manager_b, io_b, sync, "DST_JETCALO", "physical-source-b");
    the106::c0h2::beginPairAttempt(sync);
    const the106::c0h2::OccurrenceToken a =
        the106::c0h2::beginOrResumeOccurrence(
            manager_a, sync, "DST_JET", "physical-source-a", 0,
            the106::c0h2::OccurrenceTerminal::sync_skipped_behind_master,
            "synthetic source A");
    const the106::c0h2::OccurrenceToken b =
        the106::c0h2::beginOrResumeOccurrence(
            manager_b, sync, "DST_JETCALO", "physical-source-b", 0,
            the106::c0h2::OccurrenceTerminal::sync_skipped_behind_master,
            "synthetic source B");
    {
      the106::c0h2::ScopedLowLevelReadBinding read(
          io_a, a, the106::c0h2::LowLevelReadKind::full_event, 0);
      the106::c0h2::observeLowLevelRead(
          io_a, the106::c0h2::LowLevelReadOutcome::full_positive_bytes, 100);
    }
    {
      the106::c0h2::ScopedLowLevelReadBinding read(
          io_b, b, the106::c0h2::LowLevelReadKind::full_event, 0);
      the106::c0h2::observeLowLevelRead(
          io_b, the106::c0h2::LowLevelReadOutcome::full_positive_bytes, 200);
    }
    the106::c0h2::markOccurrenceReady(manager_a, a);
    the106::c0h2::markOccurrenceReady(manager_b, b);
    the106::c0h2::markPairReady(sync);

    {
      the106::c0h2::ScopedFrameworkDelivery delivery;
      assert(the106::c0h2::currentFrameworkPairToken().valid());
      the106::c0rh::ScopedEventObservation event(1, "RecoilJets_AuAu");
      assert(event.active());
      event.bindIdentity(true, 69022, true, 17);
      event.setCentrality(true, 42.5);
      const std::string event_triggers[] = {"MBDNS2"};
      event.setActiveTriggers(event_triggers, 1);
      const std::string names[] = {
          "weta", "wphi", "weta33", "wphi33", "weta35",
          "wphi53", "et1", "e11e33", "e32e35"};
      const std::string triggers[] = {"MBDNS2"};
      const char* tags[] = {"pre", "tight"};

      // C1 amends only the unavailable natural-multiplicity prerequisite.
      // Three candidates share one real writer event.  Their producer order
      // and object identities are intentionally distinct; the separate C1
      // fixture assigns equal ET to the first two and proves deterministic
      // tie resolution without changing this bounded raw-QA schema.
      for (std::uint64_t candidate_index = 0; candidate_index < 3;
           ++candidate_index)
      {
        event.noteEncounteredCandidate();

        the106::c0h2::CandidateContext context = {};
        context.pair = the106::c0h2::currentFrameworkPairToken();
        context.delivered_event_ordinal = 1;
        context.run_valid = true;
        context.run_number = 69022;
        context.event_valid = true;
        context.event_number = 17;
        context.container_key = 3;
        context.cluster_id = 44 + candidate_index;
        context.producer_encounter_ordinal = candidate_index;
        context.module_name = "RecoilJets_AuAu";
        context.canonical_view = true;
        context.view_key = "canonical";

        const double offset = 0.01 * static_cast<double>(candidate_index);
        const double values[] = {
            0.0 + offset, -0.0 + offset,
            std::numeric_limits<double>::denorm_min() + offset,
            0.25 + offset, 0.5 + offset, 0.75 + offset,
            1.0 + offset, 1.1 + offset, 1.15 + offset};
        std::uint64_t value_bits[9];
        std::uint16_t valid_mask = 0;
        for (std::size_t index = 0; index < 9; ++index)
        {
          value_bits[index] = bits(values[index]);
          if (std::isfinite(values[index])) valid_mask |= (1U << index);
        }

        the106::c0rh::RawQACandidateObservation candidate = {};
        candidate.context = context;
        candidate.values = values;
        candidate.value_bits = value_bits;
        candidate.value_count = 9;
        candidate.valid_mask = valid_mask;
        candidate.variable_names = names;
        candidate.active_triggers = triggers;
        candidate.active_trigger_count = 1;
        candidate.photon_pt_slice = 0;
        candidate.centrality_slice = 1;
        candidate.canonical_view = true;
        candidate.view_suffix = "canonical";
        candidate.tight_tag = the106::c0rh::RawQATightTag::tight;
        candidate.preselection_pass = true;
        candidate.direct_candidate_admitted = true;
        assert(the106::c0rh::emitRawQACandidateObservation(candidate));

        {
          the106::c0h2::ScopedCandidateContext score_context(context);
          const float features[] = {
              -0.0F, 1.0F + static_cast<float>(candidate_index)};
          const std::string feature_names[] = {"minus_zero", "one"};
          the106::c0h2::emitScoreObservation(
              the106::c0h2::ScoreStatus::valid,
              features, 2, feature_names, "model.root", "rawqa-test",
              true, 0.5F + 0.1F * static_cast<float>(candidate_index), 0);
        }

        for (std::size_t tag_index = 0; tag_index < 2; ++tag_index)
        {
          for (std::size_t index = 0; index < 9; ++index)
          {
            const std::string object = "h_ss_" + names[index] + "_" +
                                       tags[tag_index] +
                                       "_pT_5_10_cent_30_50";
            the106::c0rh::RawQAFillObservation fill = {};
            fill.context = context;
            fill.variable_name = names[index].c_str();
            fill.trigger_name = triggers[0].c_str();
            fill.tag_name = tags[tag_index];
            fill.view_suffix = "canonical";
            fill.directory_name = triggers[0].c_str();
            fill.object_name = object.c_str();
            fill.object_class = "TH1F";
            const std::string object_contract = triggers[0] + "/" + object;
            fill.object_contract_key = object_contract.c_str();
            fill.photon_pt_slice = 0;
            fill.centrality_slice = 1;
            fill.canonical_view = true;
            fill.tight_tag = the106::c0rh::RawQATightTag::tight;
            fill.value = values[index];
            fill.value_bits = value_bits[index];
            fill.weight = 1.0;
            fill.weight_bits = bits(1.0);
            fill.fill_ordinal = the106::c0rh::nextFillOrdinal();
            fill.value_valid = true;
            fill.filled = true;
            fill.sumw2_enabled = false;
            the106::c0rh::emitRawQAFillObservation(fill);
          }
        }
      }
      event.setOutcome(the106::c0rh::EventProcessOutcome::completed, 0,
                       "synthetic completion");
    }

    the106::c0h2::unbindInputManagerIO(manager_a, io_a);
    the106::c0h2::unbindInputManagerIO(manager_b, io_b);
    assert(writer.finish());
    assert(writer.serializationFailures() == 0);
  }

  assert(lineCount(directory + "/metadata.tsv") == 14);
  assert(lineCount(directory + "/score_observations.tsv") == 4);
  assert(lineCount(directory + "/score_features.tsv") == 7);
  assert(lineCount(directory + "/event_active_triggers.tsv") == 2);
  assert(lineCount(directory + "/rawqa_candidates.tsv") == 4);
  assert(lineCount(directory + "/rawqa_candidate_values.tsv") == 28);
  assert(lineCount(directory + "/rawqa_candidate_triggers.tsv") == 4);
  assert(lineCount(directory + "/rawqa_fill_witnesses.tsv") == 55);
  assert(lineCount(directory + "/completion.tsv") > 5);

  float plus = 0.0F;
  double minus = 0.0;
  assert(the106::c0r::binary32FromHex("0x80000000", plus));
  assert(std::signbit(plus));
  assert(the106::c0r::binary64FromHex("0x8000000000000000", minus));
  assert(std::signbit(minus));
  const std::string raw("space tab\tline\npercent%");
  assert(the106::c0r::percentDecode(the106::c0r::percentEncode(raw)) == raw);

  if (!requested_directory)
  {
    // The cache namespace is immutable: a second writer must fail before it
    // can register or process an event.
    bool overwrite_rejected = false;
    try
    {
      the106::c0r::DiagnosticWriter duplicate(directory, metadata);
    }
    catch (const std::runtime_error&)
    {
      overwrite_rejected = true;
    }
    assert(overwrite_rejected);

  const char* files[] = {
      "metadata.tsv", "source_observations.tsv", "score_observations.tsv",
      "score_features.tsv", "event_observations.tsv", "rawqa_candidates.tsv",
      "event_active_triggers.tsv",
      "rawqa_candidate_values.tsv", "rawqa_candidate_triggers.tsv",
      "rawqa_fill_witnesses.tsv", "completion.tsv"};
    for (std::size_t index = 0; index < sizeof(files) / sizeof(files[0]); ++index)
      assert(::unlink((directory + "/" + files[index]).c_str()) == 0);
    assert(::rmdir(directory.c_str()) == 0);
  }
  return 0;
}
