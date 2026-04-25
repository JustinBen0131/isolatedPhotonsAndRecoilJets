#include "RunPreselectionFailsFast.C"

int RunPreselectionFailsFast_WithStyle()
{
    gROOT->Macro("macros/sPhenixStyle.C");
    return RunPreselectionFailsFast();
}
