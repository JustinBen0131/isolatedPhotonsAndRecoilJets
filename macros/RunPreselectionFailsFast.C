#include "AnalyzeRecoilJets.cpp"

int RunPreselectionFailsFast()
{
    ARJ::DO_purityAndLeakageCHECKS_ONLY = true;
    ARJ::generateUEcomparisonSSQA = false;
    ARJ::skipToCentralityAndPtOverlaysWithSSQA = false;
    ARJ::SSoverlayPerVAR_processONLY = false;
    return AnalyzeRecoilJets();
}
