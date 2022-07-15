#include <iostream>
#include <vector>

void script_PRD_ECal_InDeltaPty_Genie() {

	gROOT->ProcessLine(".L PRD_ECalInDeltaPtySlices_Genie.cxx++");
	gROOT->ProcessLine("PRD_ECalInDeltaPtySlices_Genie()");

	gROOT->ProcessLine(".q");		

}
