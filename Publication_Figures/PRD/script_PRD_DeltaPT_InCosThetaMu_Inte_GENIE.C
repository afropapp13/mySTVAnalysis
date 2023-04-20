#include <iostream>
#include <vector>

void script_PRD_DeltaPT_InCosThetaMu_Inte_GENIE() {

	gROOT->ProcessLine(".L PRD_DeltaPTInCosThetaMu_InteBreakdown_GENIE.cxx++");
	gROOT->ProcessLine("PRD_DeltaPTInCosThetaMu_InteBreakdown_GENIE()");	

	gROOT->ProcessLine(".q");	

}
