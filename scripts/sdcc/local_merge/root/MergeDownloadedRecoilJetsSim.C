// MergeDownloadedRecoilJetsSim.C
//
// Local helper called by scripts/sftp_get_recoiljets_outputs.sh after SIM
// outputs are downloaded. It materializes the canonical merged SIM products
// that analysis macros consume from dataOutput/combinedSimOnly*.

#include "/Users/patsfan753/Desktop/ThesisAnalysis/macros/AnalyzeRecoilJets.h"

#include <TSystem.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{
using std::cerr;
using std::cout;
using std::string;
using std::vector;

const string kBase = "/Users/patsfan753/Desktop/ThesisAnalysis";

bool Exists(const string& path)
{
    return gSystem && !gSystem->AccessPathName(path.c_str());
}

void EnsureParent(const string& path)
{
    ARJ::EnsureParentDirForFile(path);
}

vector<string> ReadCfgList(const string& cfgListPath)
{
    vector<string> cfgs;
    std::ifstream in(cfgListPath);
    string line;
    while (std::getline(in, line))
    {
        const auto first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos) continue;
        const auto last = line.find_last_not_of(" \t\r\n");
        const string cfg = line.substr(first, last - first + 1);
        if (!cfg.empty()) cfgs.push_back(cfg);
    }
    return cfgs;
}

bool RequireInputs(const vector<string>& inputs)
{
    bool ok = true;
    for (const auto& in : inputs)
    {
        if (!Exists(in))
        {
            cerr << "[MERGE DOWNLOADED SIM][WARN] Missing input: " << in << "\n";
            ok = false;
        }
    }
    return ok;
}

bool MergeOne(const string& dataset, const string& cfg)
{
    if (dataset == "isSim")
    {
        const vector<string> inputs = {
            kBase + "/InputFiles/simPhotonJet/RecoilJets_photonjet5_ALL_" + cfg + ".root",
            kBase + "/InputFiles/simPhotonJet/RecoilJets_photonjet10_ALL_" + cfg + ".root",
            kBase + "/InputFiles/simPhotonJet/RecoilJets_photonjet20_ALL_" + cfg + ".root"
        };
        const string out = kBase + "/dataOutput/combinedSimOnly/" + cfg +
            "/photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root";
        if (!RequireInputs(inputs)) return false;
        return ARJ::BuildMergedSIMFile_PhotonSlices(
            inputs,
            {ARJ::kSigmaPhoton5_pb, ARJ::kSigmaPhoton10_pb, ARJ::kSigmaPhoton20_pb},
            out,
            ARJ::kDirSIM,
            {"photonJet5", "photonJet10", "photonJet20"}
        );
    }

    if (dataset == "isSimEmbedded")
    {
        const vector<string> inputs = {
            kBase + "/InputFiles/simEmbedded/RecoilJets_embeddedPhoton12_ALL_" + cfg + ".root",
            kBase + "/InputFiles/simEmbedded/RecoilJets_embeddedPhoton20_ALL_" + cfg + ".root"
        };
        const string out = kBase + "/dataOutput/combinedSimOnlyEMBEDDED/" + cfg +
            "/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
        if (!RequireInputs(inputs)) return false;
        return ARJ::BuildMergedSIMFile_PhotonSlices(
            inputs,
            {ARJ::kSigmaEmbeddedPhoton12_pb, ARJ::kSigmaEmbeddedPhoton20_pb},
            out,
            ARJ::kDirSIM,
            {"embeddedPhoton12", "embeddedPhoton20"}
        );
    }

    if (dataset == "isSimEmbeddedInclusive")
    {
        const vector<string> inputs = {
            kBase + "/InputFiles/InclusiveJetSIM_EMBEDDED/RecoilJets_embeddedJet12_ALL_" + cfg + ".root",
            kBase + "/InputFiles/InclusiveJetSIM_EMBEDDED/RecoilJets_embeddedJet20_ALL_" + cfg + ".root"
        };
        const string out = kBase + "/dataOutput/combinedSimOnlyEMBEDDED/" + cfg +
            "/embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root";
        if (!RequireInputs(inputs)) return false;
        return ARJ::BuildMergedSIMFile_PhotonSlices(
            inputs,
            {ARJ::kSigmaEmbeddedInclusiveJet12_pb, ARJ::kSigmaEmbeddedInclusiveJet20_pb},
            out,
            ARJ::kDirSIM,
            {"embeddedJet12", "embeddedJet20"}
        );
    }

    if (dataset == "isSimInclusive")
    {
        const vector<string> inputs = {
            kBase + "/InputFiles/InclusiveJetSIM/RecoilJets_jet5_ALL_" + cfg + ".root",
            kBase + "/InputFiles/InclusiveJetSIM/RecoilJets_jet8_ALL_" + cfg + ".root",
            kBase + "/InputFiles/InclusiveJetSIM/RecoilJets_jet12_ALL_" + cfg + ".root",
            kBase + "/InputFiles/InclusiveJetSIM/RecoilJets_jet20_ALL_" + cfg + ".root",
            kBase + "/InputFiles/InclusiveJetSIM/RecoilJets_jet30_ALL_" + cfg + ".root",
            kBase + "/InputFiles/InclusiveJetSIM/RecoilJets_jet40_ALL_" + cfg + ".root"
        };
        const string out = kBase + "/dataOutput/combinedSimOnly/" + cfg +
            "/inclusiveJet5to40_SIM/RecoilJets_jet5plus8plus12plus20plus30plus40_MERGED.root";
        if (!RequireInputs(inputs)) return false;
        return ARJ::BuildMergedSIMFile_PhotonSlices(
            inputs,
            {ARJ::kSigmaInclusiveJet5_pb, ARJ::kSigmaInclusiveJet8_pb,
             ARJ::kSigmaInclusiveJet12_pb, ARJ::kSigmaInclusiveJet20_pb,
             ARJ::kSigmaInclusiveJet30_pb, ARJ::kSigmaInclusiveJet40_pb},
            out,
            ARJ::kDirSIM,
            {"jet5", "jet8", "jet12", "jet20", "jet30", "jet40"}
        );
    }

    cerr << "[MERGE DOWNLOADED SIM][ERROR] Unsupported dataset: " << dataset << "\n";
    return false;
}
}

void MergeDownloadedRecoilJetsSim(const char* datasetArg = "", const char* cfgListPathArg = "")
{
    const string dataset = datasetArg ? datasetArg : "";
    const string cfgListPath = cfgListPathArg ? cfgListPathArg : "";
    const vector<string> cfgs = ReadCfgList(cfgListPath);

    if (dataset.empty() || cfgs.empty())
    {
        cerr << "[MERGE DOWNLOADED SIM][ERROR] Usage: MergeDownloadedRecoilJetsSim(\"isSim\", \"/path/to/cfgs.txt\")\n";
        gSystem->Exit(2);
        return;
    }

    cout << "[MERGE DOWNLOADED SIM] dataset=" << dataset
         << " cfgCount=" << cfgs.size() << "\n";

    int okCount = 0;
    int failCount = 0;
    for (const auto& cfg : cfgs)
    {
        cout << "\n[MERGE DOWNLOADED SIM] cfg=" << cfg << "\n";
        if (MergeOne(dataset, cfg)) ++okCount;
        else ++failCount;
    }

    cout << "\n[MERGE DOWNLOADED SIM] complete: ok=" << okCount
         << " failed=" << failCount << "\n";

    if (failCount > 0) gSystem->Exit(10);
}
