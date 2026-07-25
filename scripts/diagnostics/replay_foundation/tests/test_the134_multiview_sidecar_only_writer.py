#!/usr/bin/env python3
"""Local contract tests for THE-134 sidecar-only replay validation."""

from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


HERE = Path(__file__).resolve()
REPOSITORY = HERE.parents[4]
FOUNDATION = REPOSITORY / "src" / "RJReplayFoundationV1.h"
RUNTIME = REPOSITORY / "src" / "RJReplayRuntimeV1.h"
PP_FACADE = REPOSITORY / "src" / "RecoilJets.cc"
PP_HEADER = REPOSITORY / "src" / "RecoilJets.h"
AUAU_FACADE = REPOSITORY / "src_AuAu" / "RecoilJets_AuAu.cc"
AUAU_HEADER = REPOSITORY / "src_AuAu" / "RecoilJets_AuAu.h"
CONTROLLER = (
    REPOSITORY
    / "scripts"
    / "sdcc"
    / "workflows"
    / "diagnostics"
    / "resolve_the134_full_multiview_extraction.py"
)


def root_config_path() -> Path | None:
    candidates = [
        shutil.which("root-config"),
        str(Path(sys.executable).with_name("root-config")),
        str(REPOSITORY.parent / "analysis" / "env" / "bin" / "root-config"),
        str(REPOSITORY.parents[1] / "analysis" / "env" / "bin" / "root-config"),
        str(
            Path.home()
            / "Desktop"
            / "analysis"
            / "env"
            / "bin"
            / "root-config"
        ),
    ]
    rootsys = os.environ.get("ROOTSYS")
    if rootsys:
        candidates.append(str(Path(rootsys) / "bin" / "root-config"))
    for candidate in candidates:
        if candidate and Path(candidate).is_file() and os.access(candidate, os.X_OK):
            return Path(candidate)
    return None


def root_config(config: Path, option: str) -> str:
    result = subprocess.run(
        [str(config), option],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return result.stdout.strip()


def assert_static_contract(source: str, fragments: tuple[str, ...]) -> None:
    for fragment in fragments:
        if fragment not in source:
            raise AssertionError(f"missing static contract fragment: {fragment}")


CPP_FIXTURE = r"""
#include "RJReplayRuntimeV1.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TObject.h>

#include <cstdlib>
#include <iostream>
#include <string>

using namespace RJReplayFoundationV1;
using namespace RJReplayRuntimeV1;

namespace
{
void freezeMetadata()
{
  const std::string digest(64,'a');
  setenv("RJ_REPLAY_SCHEMA_SHA256",digest.c_str(),1);
  setenv("RJ_REPLAY_SEMANTIC_SHA256",digest.c_str(),1);
  setenv("RJ_REPLAY_SOURCE_SHA256",digest.c_str(),1);
  setenv("RJ_REPLAY_MODEL_SHA256",digest.c_str(),1);
  setenv("RJ_REPLAY_CONFIG_SHA256",digest.c_str(),1);
  setenv("RJ_REPLAY_CODE_SHA256",digest.c_str(),1);
}

SourceOccurrenceRow source()
{
  SourceOccurrenceRow row;
  row.id=makeIdentity("the134-sidecar-only-source");
  row.lane="fixture";
  row.dataset="fixture";
  row.sample="fixture";
  row.period="fixture";
  row.si_di_role="SI";
  row.ownership_state="source_role_frozen";
  return row;
}

EventBundle validBundle(const SourceOccurrenceRow& sourceRow)
{
  EventBundle bundle;
  bundle.event.id=makeIdentity("the134-sidecar-only-valid-event");
  bundle.event.source_id=sourceRow.id;
  bundle.event.event_sequence=1;

  PhotonCandidateRow candidate;
  candidate.id=makeIdentity("the134-sidecar-only-valid-candidate");
  candidate.event_id=bundle.event.id;
  bundle.candidates.push_back(candidate);

  ModelEvaluationRow model;
  model.candidate_id=candidate.id;
  model.model_id=makeIdentity("the134-sidecar-only-valid-model");
  model.shower_definition_id="H70";
  model.shower_semantic_sha256=std::string(64,'a');
  bundle.models.push_back(model);

  ShowerCellRow cell;
  cell.candidate_id=candidate.id;
  cell.tower_eta_index=10;
  cell.tower_phi_index=20;
  cell.tower_key=1;
  cell.grid_membership_bitmask=1;
  bundle.shower_cells.push_back(cell);

  ShowerFeatureViewRow view;
  view.candidate_id=candidate.id;
  view.definition_id=makeIdentity("the134-sidecar-only-valid-view");
  view.definition_name="H70";
  view.semantic_sha256=std::string(64,'a');
  view.center_eta_index=10;
  view.center_phi_index=20;
  view.raw_center_eta=0.1;
  view.raw_center_phi=0.2;
  bundle.shower_feature_views.push_back(view);

  IsolationConstituentRow isolationConstituent;
  isolationConstituent.candidate_id=candidate.id;
  isolationConstituent.constituent_id=
      makeIdentity("the134-sidecar-only-valid-isolation-constituent");
  bundle.isolation_constituents.push_back(isolationConstituent);

  IsolationWitnessRow isolationWitness;
  isolationWitness.candidate_id=candidate.id;
  isolationWitness.isolation_id=
      makeIdentity("the134-sidecar-only-valid-isolation-witness");
  bundle.isolation_witnesses.push_back(isolationWitness);

  JetRow jet;
  jet.id=makeIdentity("the134-sidecar-only-valid-jet");
  jet.event_id=bundle.event.id;
  bundle.jets.push_back(jet);

  JetConstituentRow jetConstituent;
  jetConstituent.jet_id=jet.id;
  jetConstituent.constituent_id=
      makeIdentity("the134-sidecar-only-valid-jet-constituent");
  bundle.jet_constituents.push_back(jetConstituent);

  PhotonJetPairRow pair;
  pair.id=makeIdentity("the134-sidecar-only-valid-pair");
  pair.event_id=bundle.event.id;
  pair.candidate_id=candidate.id;
  pair.jet_id=jet.id;
  bundle.pairs.push_back(pair);

  TruthPhotonRow truthPhoton;
  truthPhoton.id=makeIdentity("the134-sidecar-only-valid-truth-photon");
  truthPhoton.event_id=bundle.event.id;
  bundle.truth_photons.push_back(truthPhoton);

  TruthJetRow truthJet;
  truthJet.id=makeIdentity("the134-sidecar-only-valid-truth-jet");
  truthJet.event_id=bundle.event.id;
  bundle.truth_jets.push_back(truthJet);

  RecoTruthLinkRow link;
  link.id=makeIdentity("the134-sidecar-only-valid-link");
  link.reco_id=candidate.id;
  link.truth_id=truthPhoton.id;
  link.reco_type=static_cast<int>(RecoTruthType::PHOTON);
  link.truth_type=static_cast<int>(RecoTruthType::PHOTON);
  link.link_class=static_cast<int>(LinkClass::MATCH);
  bundle.links.push_back(link);

  WeightComponentRow weight;
  weight.target_id=candidate.id;
  weight.component_type="fixture";
  weight.application_count=1;
  bundle.weights.push_back(weight);

  EventDisplaySnapshotRow snapshot;
  snapshot.id=makeIdentity("the134-sidecar-only-valid-snapshot");
  snapshot.event_id=bundle.event.id;
  bundle.snapshots.push_back(snapshot);
  return bundle;
}

EventBundle invalidCollectionBundle(const SourceOccurrenceRow& sourceRow,
                                    int collectionIndex)
{
  EventBundle bundle=validBundle(sourceRow);
  const Identity128 missing=makeIdentity(
      "the134-sidecar-only-missing-"+std::to_string(collectionIndex));
  switch(collectionIndex)
  {
    case 0: bundle.event.source_id=missing; break;
    case 1: bundle.candidates[0].event_id=missing; break;
    case 2: bundle.models[0].candidate_id=missing; break;
    case 3: bundle.shower_cells[0].tower_eta_index=-1; break;
    case 4: bundle.shower_feature_views[0].center_phi_index=-1; break;
    case 5: bundle.isolation_constituents[0].candidate_id=missing; break;
    case 6: bundle.isolation_witnesses[0].candidate_id=missing; break;
    case 7: bundle.jets[0].event_id=missing; break;
    case 8: bundle.jet_constituents[0].jet_id=missing; break;
    case 9: bundle.pairs[0].candidate_id=missing; break;
    case 10: bundle.truth_photons[0].event_id=missing; break;
    case 11: bundle.truth_jets[0].event_id=missing; break;
    case 12: bundle.links[0].reco_type=static_cast<int>(RecoTruthType::NONE); break;
    case 13: bundle.weights[0].target_id=Identity128{}; break;
    case 14: bundle.snapshots[0].event_id=missing; break;
    default: break;
  }
  return bundle;
}

EventBundle invalidBundle(const SourceOccurrenceRow& sourceRow)
{
  EventBundle bundle;
  bundle.event.id=makeIdentity("the134-sidecar-only-invalid-event");
  bundle.event.source_id=sourceRow.id;
  bundle.event.event_sequence=2;
  PhotonCandidateRow candidate;
  candidate.id=makeIdentity("the134-sidecar-only-orphan-candidate");
  candidate.event_id=makeIdentity("the134-sidecar-only-missing-event");
  bundle.candidates.push_back(candidate);
  return bundle;
}

EventBundle duplicateBundle(const SourceOccurrenceRow& sourceRow)
{
  EventBundle bundle=validBundle(sourceRow);
  bundle.event.id=makeIdentity("the134-sidecar-only-duplicate-event");
  bundle.candidates[0].event_id=bundle.event.id;
  bundle.candidates.push_back(bundle.candidates[0]);
  return bundle;
}

bool initialize(Runtime& runtime,
                TFile& file,
                const SourceOccurrenceRow& sourceRow,
                WriterMode mode,
                std::string& error)
{
  if(mode==WriterMode::SERIALIZE)
    return runtime.initialize(&file,sourceRow,&error);
  return runtime.initialize(&file,sourceRow,mode,&error);
}

bool runValid(const std::string& path,WriterMode mode)
{
  TFile file(path.c_str(),"RECREATE");
  Runtime runtime;
  const auto sourceRow=source();
  std::string error;
  if(!initialize(runtime,file,sourceRow,mode,error))
  {
    std::cerr<<"valid initialize failed: "<<error<<std::endl;
    return false;
  }
  if(file.GetCompressionAlgorithm()!=
         static_cast<int>(ROOT::RCompressionSetting::EAlgorithm::kZSTD)||
     file.GetCompressionLevel()!=5)
  {
    std::cerr<<"writer mode changed the frozen analysis-file compression"
             <<std::endl;
    return false;
  }
  const std::string initializedDigest(64,'a');
  const std::string driftedDigest(64,'b');
  setenv("RJ_REPLAY_SCHEMA_SHA256",driftedDigest.c_str(),1);
  if(runtime.metadata().schema_sha256!=initializedDigest)
  {
    std::cerr<<"runtime metadata changed after environment drift"<<std::endl;
    return false;
  }
  freezeMetadata();
  auto bundle=validBundle(sourceRow);
  if(!runtime.write(bundle,&error))
  {
    std::cerr<<"valid transaction failed: "<<error<<std::endl;
    return false;
  }
  if(!runtime.finish(&error))
  {
    std::cerr<<"valid finish failed: "<<error<<std::endl;
    return false;
  }
  error.clear();
  if(runtime.finish(&error)||error!="writer is not active")
  {
    std::cerr<<"finished writer did not reject a second finish: "<<error<<std::endl;
    return false;
  }
  file.Close();
  return true;
}

bool runInvalid(const std::string& path,WriterMode mode)
{
  TFile file(path.c_str(),"RECREATE");
  Runtime runtime;
  const auto sourceRow=source();
  std::string error;
  if(!initialize(runtime,file,sourceRow,mode,error))
  {
    std::cerr<<"invalid initialize failed: "<<error<<std::endl;
    return false;
  }
  auto bundle=invalidBundle(sourceRow);
  if(runtime.write(bundle,&error))
  {
    std::cerr<<"orphan candidate was accepted"<<std::endl;
    return false;
  }
  if(error!="candidate foreign key or identity failure")
  {
    std::cerr<<"unexpected invalid-bundle error: "<<error<<std::endl;
    return false;
  }
  file.Close();
  return true;
}

bool runDuplicate(const std::string& path,WriterMode mode)
{
  TFile file(path.c_str(),"RECREATE");
  Runtime runtime;
  const auto sourceRow=source();
  std::string error;
  if(!initialize(runtime,file,sourceRow,mode,error))
  {
    std::cerr<<"duplicate initialize failed: "<<error<<std::endl;
    return false;
  }
  auto bundle=duplicateBundle(sourceRow);
  if(runtime.write(bundle,&error))
  {
    std::cerr<<"duplicate candidate was accepted"<<std::endl;
    return false;
  }
  if(error!="candidate foreign key or identity failure")
  {
    std::cerr<<"unexpected duplicate-bundle error: "<<error<<std::endl;
    return false;
  }
  file.Close();
  return true;
}

bool rejectEveryInvalidCollection(const std::string& path,WriterMode mode)
{
  for(int collectionIndex=0;collectionIndex<15;++collectionIndex)
  {
    const std::string indexedPath=
        path+".collection-"+std::to_string(collectionIndex);
    TFile file(indexedPath.c_str(),"RECREATE");
    Runtime runtime;
    const auto sourceRow=source();
    std::string error;
    if(!initialize(runtime,file,sourceRow,mode,error))
    {
      std::cerr<<"collection initialize failed: "<<collectionIndex<<" "
               <<error<<std::endl;
      return false;
    }
    auto bundle=invalidCollectionBundle(sourceRow,collectionIndex);
    if(runtime.write(bundle,&error))
    {
      std::cerr<<"invalid collection was accepted: "<<collectionIndex
               <<std::endl;
      return false;
    }
    if(error.empty())
    {
      std::cerr<<"invalid collection lacked an error: "<<collectionIndex
               <<std::endl;
      return false;
    }
    file.Close();
  }
  return true;
}

bool rejectInvalidMetadata(const std::string& path,WriterMode mode)
{
  TFile file(path.c_str(),"RECREATE");
  Runtime runtime;
  const auto sourceRow=source();
  std::string error;
  setenv("RJ_REPLAY_SCHEMA_SHA256","not-a-sha256",1);
  const bool accepted=initialize(runtime,file,sourceRow,mode,error);
  freezeMetadata();
  file.Close();
  if(accepted)
  {
    std::cerr<<"invalid metadata was accepted"<<std::endl;
    return false;
  }
  return error=="metadata hashes must be 64 lowercase hexadecimal characters";
}

bool inspectSerialized(const std::string& path)
{
  TFile file(path.c_str(),"READ");
  auto* directory=file.GetDirectory("ReplayFoundationV1");
  if(!directory)
  {
    std::cerr<<"serialized mode did not create ReplayFoundationV1"<<std::endl;
    return false;
  }
  const char* trees[]={
    "RJSourceOccurrenceV1","RJEventV1","RJPhotonCandidateV1",
    "RJModelEvaluationV1","RJShowerCellV1","RJShowerFeatureViewV1",
    "RJIsolationConstituentV1","RJIsolationWitnessV1","RJJetV1",
    "RJJetConstituentV1","RJPhotonJetPairV1","RJTruthPhotonV1",
    "RJTruthJetV1","RJRecoTruthLinkV1","RJWeightComponentV1",
    "RJEventDisplaySnapshotV1"};
  for(const char* tree:trees)
    if(!directory->Get(tree))
    {
      std::cerr<<"serialized mode is missing tree "<<tree<<std::endl;
      return false;
    }
  return directory->Get("rj_replay_complete")!=nullptr;
}

bool inspectValidateOnly(const std::string& path)
{
  TFile file(path.c_str(),"READ");
  if(file.GetDirectory("ReplayFoundationV1"))
  {
    std::cerr<<"validate-only mode created ReplayFoundationV1"<<std::endl;
    return false;
  }
  return file.Get("rj_replay_schema")==nullptr &&
         file.Get("rj_replay_complete")==nullptr;
}
}

int main(int argc,char** argv)
{
  if(argc!=5)return 64;
  freezeMetadata();
  if(!runValid(argv[1],WriterMode::SERIALIZE))return 1;
  if(!runValid(argv[2],WriterMode::VALIDATE_ONLY))return 2;
  if(!runInvalid(argv[3],WriterMode::SERIALIZE))return 3;
  if(!runInvalid(argv[4],WriterMode::VALIDATE_ONLY))return 4;
  if(!runDuplicate(argv[3],WriterMode::SERIALIZE))return 5;
  if(!runDuplicate(argv[4],WriterMode::VALIDATE_ONLY))return 6;
  if(!rejectInvalidMetadata(argv[3],WriterMode::SERIALIZE))return 7;
  if(!rejectInvalidMetadata(argv[4],WriterMode::VALIDATE_ONLY))return 8;
  if(!rejectEveryInvalidCollection(argv[3],WriterMode::SERIALIZE))return 9;
  if(!rejectEveryInvalidCollection(argv[4],WriterMode::VALIDATE_ONLY))return 10;
  if(!inspectSerialized(argv[1]))return 11;
  if(!inspectValidateOnly(argv[2]))return 12;
  return 0;
}
"""


class TestThe134MultiviewSidecarOnlyWriter(unittest.TestCase):
    def test_root_backed_writer_modes(self) -> None:
        config = root_config_path()
        if config is None:
            self.skipTest("root-config is unavailable")
        compiler_text = root_config(config, "--cxx")
        compiler = shlex.split(compiler_text)[0] if compiler_text else "c++"
        configured_compiler = config.parent / compiler
        if not Path(compiler).is_absolute() and configured_compiler.is_file():
            compiler = str(configured_compiler)
        if shutil.which(compiler) is None and not Path(compiler).is_file():
            self.skipTest(f"ROOT C++ compiler is unavailable: {compiler}")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "writer_mode_fixture.cc"
            executable = root / "writer_mode_fixture"
            source.write_text(textwrap.dedent(CPP_FIXTURE), encoding="utf-8")
            command = [
                compiler,
                "-std=c++17",
                "-O0",
                "-I",
                str(REPOSITORY / "src"),
                *shlex.split(root_config(config, "--cflags")),
                str(source),
                "-o",
                str(executable),
                *shlex.split(root_config(config, "--libs")),
            ]
            compiled = subprocess.run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            self.assertEqual(
                compiled.returncode,
                0,
                msg=f"compile stdout:\n{compiled.stdout}\n"
                f"compile stderr:\n{compiled.stderr}",
            )
            paths = [
                root / "serialized_valid.root",
                root / "validate_only_valid.root",
                root / "serialized_invalid.root",
                root / "validate_only_invalid.root",
            ]
            environment = dict(os.environ)
            library_directory = root_config(config, "--libdir")
            for key in ("DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"):
                current = environment.get(key)
                environment[key] = (
                    f"{library_directory}{os.pathsep}{current}"
                    if current
                    else library_directory
                )
            executed = subprocess.run(
                [str(executable), *(str(path) for path in paths)],
                env=environment,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            self.assertEqual(
                executed.returncode,
                0,
                msg=f"fixture stdout:\n{executed.stdout}\n"
                f"fixture stderr:\n{executed.stderr}",
            )

    def test_facade_flag_combinations_and_marker_order_are_frozen(self) -> None:
        facade_contracts = (
            (
                PP_FACADE.read_text(encoding="utf-8"),
                PP_HEADER.read_text(encoding="utf-8"),
                (
                    'envEnabled("RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1")',
                    "!m_replayFoundationEnabled || !multiviewTrainingEnabled",
                    "!m_ppPhotonIDExtractOnly || !m_ppPhotonIDTrainingTreeEnabled",
                    "!m_ppPhotonIDTrainingTreeEnabled ||\n       m_isAuAu",
                    "RJReplayFoundationV1::WriterMode::VALIDATE_ONLY",
                ),
            ),
            (
                AUAU_FACADE.read_text(encoding="utf-8"),
                AUAU_HEADER.read_text(encoding="utf-8"),
                (
                    'envEnabled("RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1")',
                    "!m_replayFoundationEnabled||!multiviewTrainingEnabled",
                    "!m_auauBDTExtractOnly||!m_auauBDTTrainingTreeEnabled",
                    "m_auauBDTTrainingTreeEnabled||!m_isAuAu",
                    "RJReplayFoundationV1::WriterMode::VALIDATE_ONLY",
                ),
            ),
        )
        markers = (
            '{"rj_the134_multiview_sidecar_only_v1","1"}',
            '{"rj_replay_transaction_state","CONSTRUCTED_AND_VALIDATED"}',
            '{"rj_replay_serialization_state","DISABLED"}',
            '{"rj_replay_cache_applicability","NOT_APPLICABLE"}',
            "m_replayRuntime->metadata()",
            '{"rj_replay_schema_sha256"',
            "replayMetadata.schema_sha256",
            '{"rj_replay_semantic_sha256"',
            "replayMetadata.semantic_sha256",
            '{"rj_replay_source_sha256"',
            "replayMetadata.source_sha256",
            '{"rj_replay_model_sha256"',
            "replayMetadata.model_sha256",
            '{"rj_replay_config_sha256"',
            "replayMetadata.config_sha256",
            '{"rj_replay_code_sha256"',
            "replayMetadata.code_sha256",
        )
        for index, (source, header, fragments) in enumerate(facade_contracts):
            with self.subTest(facade=index):
                assert_static_contract(source, fragments + markers)
                self.assertIn(
                    "bool m_the134MultiviewSidecarOnly = false;", header
                )
                self.assertEqual(source.count(markers[0]), 1)
                replay_finish = source.index(
                    "m_replayRuntime->finish(&replayError)"
                )
                sidecar_finish = source.index(
                    "m_photonTrainingViewRuntime->finish(&trainingError)"
                )
                first_marker = source.index(markers[0])
                self.assertLess(replay_finish, sidecar_finish)
                self.assertLess(sidecar_finish, first_marker)
                conditional = source.rfind(
                    "if (m_the134MultiviewSidecarOnly)", 0, first_marker
                )
                self.assertGreater(conditional, sidecar_finish)
                for fragment in fragments + markers:
                    with self.subTest(facade=index, mutation=fragment):
                        mutated = source.replace(fragment, "")
                        with self.assertRaises(AssertionError):
                            assert_static_contract(
                                mutated, fragments + markers
                            )

    def test_controller_sidecar_only_profile_is_frozen(self) -> None:
        source = CONTROLLER.read_text(encoding="utf-8")
        fragments = (
            '"RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1": "1"',
            '"artifact_profile": "THE134_MULTIVIEW_SIDECAR_ONLY_V1"',
            '"ANALYSIS_AND_LEGACY_TRAINING_WITH_VALIDATION_MARKERS"',
            '"replay_transaction": "CONSTRUCTED_AND_VALIDATED"',
            '"replay_serialization": "DISABLED"',
            '"cache_replay_applicability": "NOT_APPLICABLE"',
            '"rj_replay_schema_sha256": "RJ_REPLAY_SCHEMA_SHA256"',
            '"rj_replay_semantic_sha256": "RJ_REPLAY_SEMANTIC_SHA256"',
            '"rj_replay_source_sha256": "RJ_REPLAY_SOURCE_SHA256"',
            '"rj_replay_model_sha256": "RJ_REPLAY_MODEL_SHA256"',
            '"rj_replay_config_sha256": "RJ_REPLAY_CONFIG_SHA256"',
            '"rj_replay_code_sha256": "RJ_REPLAY_CODE_SHA256"',
            '"full_training_authority": 0',
            "def validate_sidecar_only_contract(",
        )
        assert_static_contract(source, fragments)
        for fragment in fragments:
            with self.subTest(mutation=fragment):
                mutated = source.replace(fragment, "")
                with self.assertRaises(AssertionError):
                    assert_static_contract(mutated, fragments)

    def test_runtime_write_body_remains_unchanged_by_mode(self) -> None:
        source = RUNTIME.read_text(encoding="utf-8")
        write_body = source[
            source.index("  bool write(EventBundle& bundle"):
            source.index(
                "  bool finish(std::string* error = nullptr)",
                source.index("  bool write(EventBundle& bundle"),
            )
        ]
        self.assertNotIn("WriterMode", write_body)
        self.assertNotIn("VALIDATE_ONLY", write_body)
        runtime_fragments = (
            "m_writer.fill(bundle.event, error)",
            "bundle.candidates) if (!m_writer.fill(row, error))",
            "bundle.models) if (!m_writer.fill(row, error))",
            "bundle.shower_cells) if (!m_writer.fill(row, error))",
            "bundle.shower_feature_views) if (!m_writer.fill(row, error))",
            "bundle.isolation_constituents) if (!m_writer.fill(row, error))",
            "bundle.isolation_witnesses) if (!m_writer.fill(row, error))",
            "bundle.jets) if (!m_writer.fill(row, error))",
            "bundle.jet_constituents) if (!m_writer.fill(row, error))",
            "bundle.pairs) if (!m_writer.fill(row, error))",
            "bundle.truth_photons) if (!m_writer.fill(row, error))",
            "bundle.truth_jets) if (!m_writer.fill(row, error))",
            "bundle.links) if (!m_writer.fill(row, error))",
            "bundle.weights) if (!m_writer.fill(row, error))",
            "bundle.snapshots) if (!m_writer.fill(row, error))",
        )
        assert_static_contract(write_body, runtime_fragments)
        for fragment in runtime_fragments:
            with self.subTest(runtime_mutation=fragment):
                mutated = write_body.replace(fragment, "")
                with self.assertRaises(AssertionError):
                    assert_static_contract(mutated, runtime_fragments)
        foundation = FOUNDATION.read_text(encoding="utf-8")
        self.assertEqual(foundation.count("if(validationOnly())return true;"), 16)
        self.assertIn("enum class WriterMode", FOUNDATION.read_text())


if __name__ == "__main__":
    unittest.main()
