#include "AnalyzeRecoilJets.h"

namespace ARJ {
#include "AnalyzeRecoilJets_RunIsoQA_UEComparisons_AuAu.cpp"
}

int RunSSQAFast()
{
    ARJ::generateUEcomparisonSSQA = true;
    ARJ::skipToCentralityAndPtOverlaysWithSSQA = false;
    ARJ::SSoverlayPerVAR_processONLY = true;
    ARJ::RunIsoQA_UEComparisons_AuAu(0);
    return 0;
}
