#include <iostream>
#include <vector>

void script_PRD_DeltaPT_InDeltaAlphaT_Genie() {

	gROOT->ProcessLine(".L PRD_DeltaPTInDeltaAlphaTSlices_Genie.cxx++");
	gROOT->ProcessLine("PRD_DeltaPTInDeltaAlphaTSlices_Genie()");	

	gROOT->ProcessLine(".q");	

}
