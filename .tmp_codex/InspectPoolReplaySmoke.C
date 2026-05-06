#include <TFile.h>
#include <TH1.h>
#include <TKey.h>

#include <iostream>

void InspectPoolReplaySmoke(const char* path = ".tmp_codex/pool_smoke/replay.root")
{
  TFile f(path, "READ");
  TH1* a = dynamic_cast<TH1*>(f.Get("SIM/h_isIsolated_isTight_pT_20_22"));
  TH1* c = dynamic_cast<TH1*>(f.Get("SIM/h_isIsolated_notTight_pT_20_22"));
  TH1* sig = dynamic_cast<TH1*>(f.Get("SIM/h_sigABCD_MC_pT_20_22"));
  TH1* cnt = dynamic_cast<TH1*>(f.Get("SIM/cnt_SIM"));
  std::cout << "cnt_SIM=" << (cnt ? cnt->GetBinContent(1) : -1) << "\n";
  std::cout << "A_20_22=" << (a ? a->GetBinContent(1) : -1) << "\n";
  std::cout << "C_20_22=" << (c ? c->GetBinContent(1) : -1) << "\n";
  std::cout << "sigABCD_C=" << (sig ? sig->GetBinContent(3) : -1) << "\n";
  if (!a)
  {
    auto* dir = f.GetDirectory("SIM");
    if (dir)
    {
      TIter next(dir->GetListOfKeys());
      while (TKey* key = static_cast<TKey*>(next()))
      {
        std::string name = key->GetName();
        if (name.find("isIsolated") != std::string::npos)
          std::cout << "ABCD key: " << name << "\n";
      }
    }
  }
}
