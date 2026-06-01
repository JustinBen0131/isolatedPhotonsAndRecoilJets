// File: ExtractDoNotScaleHistograms.C
//
// Usage:
//   root -l -q 'ExtractDoNotScaleHistograms.C()'
//
// What this does:
//   - Opens the source ROOT file
//   - Scans all top-level directories
//   - Copies only histograms named
//       h_maxEnergyClus_NewTriggerFilling_doNotScale_<directoryName>
//     into a new ROOT file
//   - Preserves the original directory names exactly, so downstream
//     turn-on scripts can pair:
//       <probeKey>  with  baseline_<probeKey>
//
// Output:
//   /Users/patsfan753/Desktop/RecoilJets_pp24_doNotScaleOnly.root
//
// Notes for whoever receives the output file:
//   - Ignore directories starting with "baseline_"
//   - For each probe directory P, use baseline_P as the denominator
//   - Ratio = probe / baseline
//   - Overlay ratios sharing the same common-basis label
//

#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TObjString.h>

#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
  std::string StripBaselinePrefix(std::string key)
  {
    const std::string prefix = "baseline_";
    if (key.rfind(prefix, 0) == 0)
    {
      key.erase(0, prefix.size());
    }
    return key;
  }

  std::string InferCommonBasisGroup(const std::string& rawKey)
  {
    const std::string key = StripBaselinePrefix(rawKey);

    if (key.find("_MBD_NS_geq_1_vtx_lt_10") != std::string::npos)
    {
      return "commonBasis_MBD_NandS_geq_1_vtx_lt_10";
    }
    if (key.find("_MBD_NS_geq_2_vtx_lt_150") != std::string::npos)
    {
      return "commonBasis_MBD_NS_geq_2_vtx_lt_150";
    }
    if (key.find("_MBD_NS_geq_2_vtx_lt_10") != std::string::npos)
    {
      return "commonBasis_MBD_NS_geq_2_vtx_lt_10";
    }
    if (key.find("_MBD_NS_geq_2") != std::string::npos)
    {
      return "commonBasis_MBD_NS_geq_2";
    }
    if (key.find("_MBD_NS_geq_1") != std::string::npos)
    {
      return "commonBasis_MBD_NandS_geq_1";
    }

    return "commonBasis_UNKNOWN";
  }

  bool CopyOneDoNotScaleHist(TFile* inFile,
                             TFile* outFile,
                             const std::string& dirName,
                             const std::string& histPrefix,
                             std::vector<std::string>& copiedDirs)
  {
    if (!inFile || !outFile || dirName.empty())
    {
      return false;
    }

    TDirectory* inDir = inFile->GetDirectory(dirName.c_str());
    if (!inDir)
    {
      return false;
    }

    const std::string histName = histPrefix + dirName;
    TH1* hIn = dynamic_cast<TH1*>(inDir->Get(histName.c_str()));
    if (!hIn)
    {
      return false;
    }

    TDirectory* outDir = outFile->GetDirectory(dirName.c_str());
    if (!outDir)
    {
      outDir = outFile->mkdir(dirName.c_str());
    }
    if (!outDir)
    {
      return false;
    }

    outDir->cd();

    TH1* hOut = dynamic_cast<TH1*>(hIn->Clone(histName.c_str()));
    if (!hOut)
    {
      outFile->cd();
      return false;
    }

    hOut->SetDirectory(outDir);
    hOut->Write("", TObject::kOverwrite);
    delete hOut;

    copiedDirs.push_back(dirName);
    outFile->cd();
    return true;
  }
}

void ExtractDoNotScaleHistograms(
    const char* inputPath =
      "/Users/patsfan753/Desktop/ThesisAnalysis/InputFiles/pp24/"
      "RecoilJets_pp_ALL_jetMinPt5_7pi_8_vz30_isoR30_isSliding.root",
    const char* outputPath =
      "/Users/patsfan753/Desktop/RecoilJets_pp24_doNotScaleOnly.root")
{
  const std::string histPrefix = "h_maxEnergyClus_NewTriggerFilling_doNotScale_";

  TFile* inFile = TFile::Open(inputPath, "READ");
  if (!inFile || inFile->IsZombie())
  {
    std::cerr << "[ERROR] Could not open input file:\n  " << inputPath << std::endl;
    if (inFile)
    {
      inFile->Close();
      delete inFile;
    }
    return;
  }

  TFile* outFile = TFile::Open(outputPath, "RECREATE");
  if (!outFile || outFile->IsZombie())
  {
    std::cerr << "[ERROR] Could not create output file:\n  " << outputPath << std::endl;
    if (outFile)
    {
      outFile->Close();
      delete outFile;
    }
    inFile->Close();
    delete inFile;
    return;
  }

  std::vector<std::string> copiedDirs;
  copiedDirs.reserve(64);

  TIter nextKey(inFile->GetListOfKeys());
  while (TObject* obj = nextKey())
  {
    TKey* key = dynamic_cast<TKey*>(obj);
    if (!key)
    {
      continue;
    }

    const std::string className = key->GetClassName();
    if (className != "TDirectoryFile" && className != "TDirectory")
    {
      continue;
    }

    const std::string dirName = key->GetName();
    CopyOneDoNotScaleHist(inFile, outFile, dirName, histPrefix, copiedDirs);
  }

  std::sort(copiedDirs.begin(), copiedDirs.end());

  std::set<std::string> groupsPresent;
  std::ostringstream readme;
  readme << "This file contains only the doNotScale histograms extracted from:\n"
         << inputPath << "\n\n"
         << "Schema preserved exactly from the source file:\n"
         << "  directory name = <probeKey> or baseline_<probeKey>\n"
         << "  histogram name = h_maxEnergyClus_NewTriggerFilling_doNotScale_<directoryName>\n\n"
         << "How to rebuild the trigger-efficiency turn-ons:\n"
         << "  1. Ignore directories whose names start with baseline_\n"
         << "  2. For each probe directory P, use baseline_P as the denominator\n"
         << "  3. Compute ratio = probe / baseline\n"
         << "  4. Overlay ratios that share the same inferred common-basis label\n\n"
         << "Common-basis inference used here:\n"
         << "  *_MBD_NS_geq_1            -> commonBasis_MBD_NandS_geq_1\n"
         << "  *_MBD_NS_geq_1_vtx_lt_10  -> commonBasis_MBD_NandS_geq_1_vtx_lt_10\n"
         << "  *_MBD_NS_geq_2            -> commonBasis_MBD_NS_geq_2\n"
         << "  *_MBD_NS_geq_2_vtx_lt_10  -> commonBasis_MBD_NS_geq_2_vtx_lt_10\n"
         << "  *_MBD_NS_geq_2_vtx_lt_150 -> commonBasis_MBD_NS_geq_2_vtx_lt_150\n\n"
         << "Copied directories (" << copiedDirs.size() << "):\n";

  for (const auto& dirName : copiedDirs)
  {
    const std::string group = InferCommonBasisGroup(dirName);
    groupsPresent.insert(group);
    readme << "  " << dirName << "    [" << group << "]\n";
  }

  readme << "\nGroups present:\n";
  for (const auto& group : groupsPresent)
  {
    readme << "  " << group << "\n";
  }

  outFile->cd();
  TObjString meta(readme.str().c_str());
  meta.Write("README_doNotScale", TObject::kOverwrite);

  outFile->Write();
  outFile->Close();
  inFile->Close();

  std::cout << "\n[OK] Wrote " << copiedDirs.size()
            << " doNotScale histogram directories to:\n  "
            << outputPath << "\n" << std::endl;

  std::cout << "Copied objects:\n";
  for (const auto& dirName : copiedDirs)
  {
    std::cout << "  " << dirName
              << "   ->   "
              << histPrefix << dirName
              << "   [" << InferCommonBasisGroup(dirName) << "]"
              << std::endl;
  }

  delete outFile;
  delete inFile;
}
