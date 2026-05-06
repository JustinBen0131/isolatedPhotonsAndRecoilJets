#include <TFile.h>
#include <TTree.h>

#include <string>
#include <vector>

void MakePoolReplaySmoke(const char* outPath = ".tmp_codex/pool_smoke/pool.root")
{
  TFile f(outPath, "RECREATE");

  int schema = 1, run = 1, evt = 1, isSim = 1, isAuAu = 0, centBin = -1, centIdx = -1;
  Long64_t eventKey = 1;
  float vz = 0.0f, weight = 1.0f;
  std::vector<std::string>* triggers = new std::vector<std::string>{"SIM"};

  TTree ev("AnalysisEventPool", "smoke event pool");
  ev.Branch("schema", &schema, "schema/I");
  ev.Branch("run", &run, "run/I");
  ev.Branch("evt", &evt, "evt/I");
  ev.Branch("eventKey", &eventKey, "eventKey/L");
  ev.Branch("isSim", &isSim, "isSim/I");
  ev.Branch("isAuAu", &isAuAu, "isAuAu/I");
  ev.Branch("centBin", &centBin, "centBin/I");
  ev.Branch("centIdx", &centIdx, "centIdx/I");
  ev.Branch("vz", &vz, "vz/F");
  ev.Branch("weight", &weight, "weight/F");
  ev.Branch("triggers", &triggers);
  ev.Fill();

  int index = 0, truthSignal = 1, truthTrackId = 42;
  float pt = 20.0f, eta = 0.1f, phi = 0.0f, energy = 21.0f, eiso = 0.5f;
  float weta = 0.45f, wphi = 0.45f, weta33 = 0.45f, wphi33 = 0.45f;
  float weta35 = 0.45f, wphi53 = 0.45f, et1 = 0.8f, et2 = 0.2f, et3 = 0.1f, et4 = 0.1f;
  float e11e33 = 0.96f, e32e35 = 0.9f;
  float ratio = 0.9f, w32 = 0.4f, w52 = 0.4f, w72 = 0.4f, meanTime = 0.0f;
  float npb = 0.8f, tightBdt = 0.9f, auauNpb = 0.8f, auauTightBdt = 0.9f;
  float truthPt = 20.0f, truthEta = 0.1f, truthPhi = 0.0f, truthIso = 0.2f;

  TTree ph("AnalysisPhotonPool", "smoke photon pool");
  ph.Branch("schema", &schema, "schema/I");
  ph.Branch("run", &run, "run/I");
  ph.Branch("evt", &evt, "evt/I");
  ph.Branch("eventKey", &eventKey, "eventKey/L");
  ph.Branch("centBin", &centBin, "centBin/I");
  ph.Branch("centIdx", &centIdx, "centIdx/I");
  ph.Branch("index", &index, "index/I");
  ph.Branch("pt", &pt, "pt/F");
  ph.Branch("eta", &eta, "eta/F");
  ph.Branch("phi", &phi, "phi/F");
  ph.Branch("energy", &energy, "energy/F");
  ph.Branch("eiso", &eiso, "eiso/F");
  ph.Branch("weta_cogx", &weta, "weta_cogx/F");
  ph.Branch("wphi_cogx", &wphi, "wphi_cogx/F");
  ph.Branch("weta33_cogx", &weta33, "weta33_cogx/F");
  ph.Branch("wphi33_cogx", &wphi33, "wphi33_cogx/F");
  ph.Branch("weta35_cogx", &weta35, "weta35_cogx/F");
  ph.Branch("wphi53_cogx", &wphi53, "wphi53_cogx/F");
  ph.Branch("et1", &et1, "et1/F");
  ph.Branch("et2", &et2, "et2/F");
  ph.Branch("et3", &et3, "et3/F");
  ph.Branch("et4", &et4, "et4/F");
  ph.Branch("e11_over_e33", &e11e33, "e11_over_e33/F");
  ph.Branch("e32_over_e35", &e32e35, "e32_over_e35/F");
  ph.Branch("e11_over_e22", &ratio, "e11_over_e22/F");
  ph.Branch("e11_over_e13", &ratio, "e11_over_e13/F");
  ph.Branch("e11_over_e15", &ratio, "e11_over_e15/F");
  ph.Branch("e11_over_e17", &ratio, "e11_over_e17/F");
  ph.Branch("e11_over_e31", &ratio, "e11_over_e31/F");
  ph.Branch("e11_over_e51", &ratio, "e11_over_e51/F");
  ph.Branch("e11_over_e71", &ratio, "e11_over_e71/F");
  ph.Branch("e22_over_e33", &ratio, "e22_over_e33/F");
  ph.Branch("e22_over_e35", &ratio, "e22_over_e35/F");
  ph.Branch("e22_over_e37", &ratio, "e22_over_e37/F");
  ph.Branch("e22_over_e53", &ratio, "e22_over_e53/F");
  ph.Branch("w32", &w32, "w32/F");
  ph.Branch("w52", &w52, "w52/F");
  ph.Branch("w72", &w72, "w72/F");
  ph.Branch("mean_time", &meanTime, "mean_time/F");
  ph.Branch("npb_score", &npb, "npb_score/F");
  ph.Branch("tight_bdt_score", &tightBdt, "tight_bdt_score/F");
  ph.Branch("auau_npb_score", &auauNpb, "auau_npb_score/F");
  ph.Branch("auau_tight_bdt_score", &auauTightBdt, "auau_tight_bdt_score/F");
  ph.Branch("truthSignal", &truthSignal, "truthSignal/I");
  ph.Branch("truthTrackId", &truthTrackId, "truthTrackId/I");
  ph.Branch("truthPt", &truthPt, "truthPt/F");
  ph.Branch("truthEta", &truthEta, "truthEta/F");
  ph.Branch("truthPhi", &truthPhi, "truthPhi/F");
  ph.Branch("truthIso", &truthIso, "truthIso/F");
  ph.Fill();

  std::string rKey = "r02";
  int isTruth = 0, jetIndex = 0;
  float jetPt = 12.0f, jetEta = 0.2f, jetPhi = 3.14f;
  TTree jt("AnalysisJetPool", "smoke jet pool");
  jt.Branch("schema", &schema, "schema/I");
  jt.Branch("run", &run, "run/I");
  jt.Branch("evt", &evt, "evt/I");
  jt.Branch("eventKey", &eventKey, "eventKey/L");
  jt.Branch("rKey", &rKey);
  jt.Branch("isTruth", &isTruth, "isTruth/I");
  jt.Branch("index", &jetIndex, "index/I");
  jt.Branch("pt", &jetPt, "pt/F");
  jt.Branch("eta", &jetEta, "eta/F");
  jt.Branch("phi", &jetPhi, "phi/F");
  jt.Fill();
  isTruth = 1;
  jt.Fill();

  int truthIndex = 0, truthBarcode = 7;
  TTree tr("AnalysisTruthPhotonPool", "smoke truth photon pool");
  tr.Branch("schema", &schema, "schema/I");
  tr.Branch("run", &run, "run/I");
  tr.Branch("evt", &evt, "evt/I");
  tr.Branch("eventKey", &eventKey, "eventKey/L");
  tr.Branch("index", &truthIndex, "index/I");
  tr.Branch("trackId", &truthTrackId, "trackId/I");
  tr.Branch("barcode", &truthBarcode, "barcode/I");
  tr.Branch("pt", &truthPt, "pt/F");
  tr.Branch("eta", &truthEta, "truthEta/F");
  tr.Branch("phi", &truthPhi, "truthPhi/F");
  tr.Branch("iso", &truthIso, "iso/F");
  tr.Fill();

  f.Write();
  delete triggers;
}
