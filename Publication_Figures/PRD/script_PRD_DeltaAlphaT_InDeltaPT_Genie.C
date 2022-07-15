#include <iostream>
#include <vector>

void script_PRD_DeltaAlphaT_InDeltaPT_Genie() {

	gROOT->ProcessLine(".L PRD_DeltaAlphaTInDeltaPTSlices_Genie.cxx++");
	gROOT->ProcessLine("PRD_DeltaAlphaTInDeltaPTSlices_Genie()");

	gROOT->ProcessLine(".q");		

}
