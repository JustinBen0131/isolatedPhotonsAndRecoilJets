#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd -P)"
macro="${repo_root}/macros/Fun4All_recoilJets_unified_impl.C"
tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/rj-fun4all-terminal.XXXXXX")"
trap 'rm -rf "$tmpdir"' EXIT

grep -Fq '#include <fun4all/Fun4AllSyncManager.h>' "$macro"
grep -Fq 'Fun4AllInputManager* permittedRepeatingPedestalInputManager = nullptr;' "$macro"
grep -Fq 'permittedRepeatingPedestalInputManager = pedIn;' "$macro"
grep -Fq 'RECOILJETS_FUN4ALL_STATUS_V2' "$macro"
grep -Fq 'RECOILJETS_FUN4ALL_INPUT_STATUS_V1' "$macro"
[[ "$(grep -Fc -- '->Repeat();' "$macro")" == 1 ]]

sed -n \
  '/BEGIN RJ_FUN4ALL_TERMINAL_CLASSIFIER_V1_PURE/,/END RJ_FUN4ALL_TERMINAL_CLASSIFIER_V1_PURE/p' \
  "$macro" |
  sed '/RJ_FUN4ALL_TERMINAL_CLASSIFIER_V1_PURE/d' \
  > "${tmpdir}/classifier.inc"
[[ -s "${tmpdir}/classifier.inc" ]]
sed -n \
  '/^  inline void enforce_fun4all_status(/,/^  \/\/\/ Trim whitespace/p' \
  "$macro" |
  sed '$d' \
  > "${tmpdir}/adapter.inc"
[[ -s "${tmpdir}/adapter.inc" ]]

cat > "${tmpdir}/classifier.cpp" <<'CPP'
#include <cstdlib>
#include <iostream>

namespace detail
{
CPP
cat "${tmpdir}/classifier.inc" >> "${tmpdir}/classifier.cpp"
cat >> "${tmpdir}/classifier.cpp" <<'CPP'
}

namespace
{
using detail::Fun4AllTerminalStatusClassV1;
using detail::Fun4AllTerminalStatusWitnessV1;

void require_class(
    const Fun4AllTerminalStatusWitnessV1& witness,
    const Fun4AllTerminalStatusClassV1 expected,
    const char* label)
{
  const auto observed = detail::classify_fun4all_terminal_status(witness);
  if (observed.classification != expected)
  {
    std::cerr << label << " classified as " << observed.status
              << "/" << observed.reason << std::endl;
    std::exit(1);
  }
}

Fun4AllTerminalStatusWitnessV1 verified_eof(
    const int ordinaryManagers,
    const bool repeatingPedestal,
    const int runNodeManagers)
{
  Fun4AllTerminalStatusWitnessV1 witness;
  witness.nEventsRequested = 0;
  witness.runRc = -ordinaryManagers;
  witness.endRc = 0;
  witness.eventOk = 0;
  witness.registeredInputManagers =
    ordinaryManagers + (repeatingPedestal ? 1 : 0) + runNodeManagers;
  witness.ordinaryInputManagers = ordinaryManagers;
  witness.exhaustedOrdinaryInputManagers = ordinaryManagers;
  witness.permittedRepeatingManagerExpected =
    repeatingPedestal ? 1 : 0;
  witness.permittedRepeatingManagerMatches =
    repeatingPedestal ? 1 : 0;
  witness.permittedRunNodeInputManagers = runNodeManagers;
  return witness;
}
}

int main()
{
  using Class = Fun4AllTerminalStatusClassV1;

  Fun4AllTerminalStatusWitnessV1 eventOk;
  require_class(eventOk, Class::kEventOk, "unchanged EVENT_OK path");

  const auto pp = verified_eof(5, true, 1);
  require_class(
      pp,
      Class::kVerifiedMultiInputEof,
      "pp five-manager EOF with run-node provider");

  const auto auau = verified_eof(4, false, 1);
  require_class(
      auau,
      Class::kVerifiedMultiInputEof,
      "AuAu four-manager EOF with run-node provider");

  auto mutation = pp;
  mutation.nEventsRequested = 1;
  require_class(mutation, Class::kFail, "bounded negative run");

  mutation = pp;
  mutation.endRc = -1;
  require_class(mutation, Class::kFail, "nonzero End");

  mutation = pp;
  mutation.abortProcessingCount = 1;
  require_class(mutation, Class::kFail, "ABORTPROCESSING statistic");

  mutation = pp;
  mutation.abortRunCount = 1;
  require_class(mutation, Class::kFail, "ABORTRUN statistic");

  mutation = pp;
  mutation.runRc = -4;
  require_class(mutation, Class::kFail, "EOF sum mismatch");

  mutation = pp;
  mutation.exhaustedOrdinaryInputManagers = 4;
  require_class(mutation, Class::kFail, "incomplete exhaustion");

  mutation = pp;
  mutation.openOrdinaryInputManagers = 1;
  require_class(mutation, Class::kFail, "ordinary manager still open");

  mutation = pp;
  mutation.nonemptyOrdinaryFileLists = 1;
  require_class(mutation, Class::kFail, "ordinary file list not empty");

  mutation = pp;
  mutation.permittedRepeatingManagerMatches = 0;
  require_class(mutation, Class::kFail, "missing repeating pointer");

  mutation = pp;
  mutation.permittedRepeatingManagerMatches = 2;
  require_class(mutation, Class::kFail, "duplicate repeating pointer");

  mutation = pp;
  mutation.registeredInputManagers = 8;
  require_class(mutation, Class::kFail, "unaccounted input manager");

  mutation = pp;
  mutation.permittedRunNodeInputManagers = -1;
  require_class(mutation, Class::kFail, "negative run-node accounting");

  mutation = auau;
  mutation.ordinaryInputManagers = 0;
  mutation.exhaustedOrdinaryInputManagers = 0;
  mutation.registeredInputManagers = 0;
  mutation.runRc = 0;
  mutation.endRc = -1;
  require_class(mutation, Class::kFail, "zero-manager false EOF");

  mutation = auau;
  mutation.openOrdinaryInputManagers = 1;
  require_class(
      mutation,
      Class::kFail,
      "unknown negative with unexhausted manager");

  std::cout << "FUN4ALL_TERMINAL_CLASSIFIER_V1_TEST_PASS"
            << " positives=3 negative_mutations=14" << std::endl;
  return 0;
}
CPP

"${CXX:-c++}" \
  -std=c++17 -Wall -Wextra -Werror \
  "${tmpdir}/classifier.cpp" \
  -o "${tmpdir}/classifier_test"
"${tmpdir}/classifier_test"

cat > "${tmpdir}/adapter.cpp" <<'CPP'
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace Fun4AllReturnCodes
{
constexpr int ABORTPROCESSING = -4;
constexpr int ABORTRUN = -2;
constexpr int EVENT_OK = 0;
}

class Fun4AllInputManager
{
 public:
  Fun4AllInputManager(
      const bool open,
      const bool empty,
      const std::string& name = "FIXTURE_INPUT")
    : m_open(open)
    , m_empty(empty)
    , m_name(name)
  {
  }
  virtual ~Fun4AllInputManager() = default;
  bool IsOpen() const { return m_open; }
  bool FileListEmpty() const { return m_empty; }
  const std::string& Name() const { return m_name; }

 private:
  bool m_open = false;
  bool m_empty = true;
  std::string m_name;
};

class Fun4AllRunNodeInputManager : public Fun4AllInputManager
{
 public:
  using Fun4AllInputManager::Fun4AllInputManager;
};

class Fun4AllSyncManager
{
 public:
  const std::vector<Fun4AllInputManager*>& GetInputManagers() const
  {
    return managers;
  }
  std::vector<Fun4AllInputManager*> managers;
};

class Fun4AllServer
{
 public:
  int retcodestats(const int code)
  {
    if (code == Fun4AllReturnCodes::ABORTPROCESSING)
    {
      return abortProcessing;
    }
    if (code == Fun4AllReturnCodes::ABORTRUN)
    {
      return abortRun;
    }
    return 0;
  }
  Fun4AllSyncManager* getSyncManager() { return sync; }

  Fun4AllSyncManager* sync = nullptr;
  int abortProcessing = 0;
  int abortRun = 0;
};

class MockSystem
{
 public:
  void Exit(const int value) { exitCode = value; }
  int exitCode = 0;
};

MockSystem mockSystem;
MockSystem* gSystem = &mockSystem;

namespace detail
{
CPP
cat "${tmpdir}/classifier.inc" >> "${tmpdir}/adapter.cpp"
cat "${tmpdir}/adapter.inc" >> "${tmpdir}/adapter.cpp"
cat >> "${tmpdir}/adapter.cpp" <<'CPP'
}

namespace
{
void require_pass(
    const char* label,
    Fun4AllServer* server,
    const int nEvents,
    const int runRc,
    const int endRc,
    const Fun4AllInputManager* repeating)
{
  try
  {
    detail::enforce_fun4all_status(
        label, server, nEvents, runRc, endRc, repeating);
  }
  catch (const std::exception& error)
  {
    std::cerr << label << " unexpectedly failed: " << error.what()
              << std::endl;
    std::exit(1);
  }
}

void require_fail(
    const char* label,
    Fun4AllServer* server,
    const int nEvents,
    const int runRc,
    const int endRc,
    const Fun4AllInputManager* repeating)
{
  mockSystem.exitCode = 0;
  try
  {
    detail::enforce_fun4all_status(
        label, server, nEvents, runRc, endRc, repeating);
  }
  catch (const std::exception&)
  {
    if (mockSystem.exitCode != 90)
    {
      std::cerr << label << " did not propagate exit 90" << std::endl;
      std::exit(1);
    }
    return;
  }
  std::cerr << label << " unexpectedly passed" << std::endl;
  std::exit(1);
}
}

int main()
{
  require_pass("event-ok", nullptr, 1, 0, 0, nullptr);

  Fun4AllInputManager exhausted(false, true);
  Fun4AllInputManager pedestal(true, false);
  Fun4AllRunNodeInputManager runNode(true, false, "DST_GEO");
  Fun4AllSyncManager ppSync;
  ppSync.managers = {
      &exhausted, &exhausted, &exhausted, &exhausted, &exhausted,
      &pedestal, &runNode};
  Fun4AllServer ppServer;
  ppServer.sync = &ppSync;
  require_pass("pp-eof", &ppServer, 0, -5, 0, &pedestal);

  Fun4AllSyncManager auauSync;
  auauSync.managers = {
      &exhausted, &exhausted, &exhausted, &exhausted, &runNode};
  Fun4AllServer auauServer;
  auauServer.sync = &auauSync;
  require_pass("auau-eof", &auauServer, 0, -4, 0, nullptr);

  Fun4AllInputManager substitute(true, false);
  require_fail(
      "pointer-substitution", &ppServer, 0, -5, 0, &substitute);

  ppServer.abortProcessing = 1;
  require_fail(
      "abort-statistic", &ppServer, 0, -5, 0, &pedestal);

  Fun4AllInputManager unknownOpen(true, false, "UNKNOWN_OPEN_INPUT");
  Fun4AllSyncManager unknownSync;
  unknownSync.managers = {
      &exhausted, &exhausted, &exhausted, &exhausted, &unknownOpen};
  Fun4AllServer unknownServer;
  unknownServer.sync = &unknownSync;
  require_fail(
      "unknown-open-ordinary", &unknownServer, 0, -4, 0, nullptr);

  std::cout << "FUN4ALL_TERMINAL_ADAPTER_V1_TEST_PASS"
            << " positives=3 negative_mutations=3" << std::endl;
  return 0;
}
CPP

"${CXX:-c++}" \
  -std=c++17 -Wall -Wextra -Werror \
  "${tmpdir}/adapter.cpp" \
  -o "${tmpdir}/adapter_test"
"${tmpdir}/adapter_test" 2>/dev/null
