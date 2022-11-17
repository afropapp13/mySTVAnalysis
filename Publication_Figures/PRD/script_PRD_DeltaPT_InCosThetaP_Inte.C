#include <iostream>
#include <vector>

void script_PRD_DeltaPT_InCosThetaP_Inte() {

	gROOT->ProcessLine(".L PRD_DeltaPTInCosThetaP_InteBreakdown.cxx++");
	gROOT->ProcessLine("PRD_DeltaPTInCosThetaP_InteBreakdown()");	

	gROOT->ProcessLine(".q");	

}
