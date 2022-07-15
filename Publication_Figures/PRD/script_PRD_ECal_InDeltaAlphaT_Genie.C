#include <iostream>
#include <vector>

void script_PRD_ECal_InDeltaAlphaT_Genie() {

	gROOT->ProcessLine(".L PRD_ECalInDeltaAlphaTSlices_Genie.cxx++");
	gROOT->ProcessLine("PRD_ECalInDeltaAlphaTSlices_Genie()");

	gROOT->ProcessLine(".q");		

}
