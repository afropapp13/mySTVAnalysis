#include <iostream>
#include <vector>

void script_PRD_DeltaAlphaT_InCosThetaP_Inte() {

	gROOT->ProcessLine(".L PRD_DeltaAlphaTInCosThetaP_InteBreakdown.cxx++");
	gROOT->ProcessLine("PRD_DeltaAlphaTInCosThetaP_InteBreakdown()");	

	gROOT->ProcessLine(".q");	

}
