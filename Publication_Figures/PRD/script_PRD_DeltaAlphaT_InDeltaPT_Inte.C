#include <iostream>
#include <vector>

void script_PRD_DeltaAlphaT_InDeltaPT_Inte() {

	gROOT->ProcessLine(".L PRD_DeltaAlphaTInDeltaPT_InteBreakdown.cxx++");
	gROOT->ProcessLine("PRD_DeltaAlphaTInDeltaPT_InteBreakdown()");

	//gROOT->ProcessLine(".q");		

}
