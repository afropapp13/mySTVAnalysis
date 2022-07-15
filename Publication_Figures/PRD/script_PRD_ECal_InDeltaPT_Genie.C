#include <iostream>
#include <vector>

void script_PRD_ECal_InDeltaPT_Genie() {

	gROOT->ProcessLine(".L PRD_ECalInDeltaPTSlices_Genie.cxx++");
	gROOT->ProcessLine("PRD_ECalInDeltaPTSlices_Genie()");

	gROOT->ProcessLine(".q");		

}
