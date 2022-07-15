#include <iostream>
#include <vector>

void script_PRD_ECal_InDeltaPT_Gene() {

	gROOT->ProcessLine(".L PRD_ECalInDeltaPTSlices_Gene.cxx++");
	gROOT->ProcessLine("PRD_ECalInDeltaPTSlices_Gene()");

	gROOT->ProcessLine(".q");		

}
