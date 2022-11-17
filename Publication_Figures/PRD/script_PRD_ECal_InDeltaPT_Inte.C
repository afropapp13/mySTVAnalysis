#include <iostream>
#include <vector>

void script_PRD_ECal_InDeltaPT_Inte() {

	gROOT->ProcessLine(".L PRD_ECalInDeltaPT_InteBreakdown.cxx++");
	gROOT->ProcessLine("PRD_ECalInDeltaPT_InteBreakdown()");

	//gROOT->ProcessLine(".q");		

}
