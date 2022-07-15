#include <iostream>
#include <vector>

void script_PRD_DeltaAlphaT_InDeltaPT_FSI() {

	gROOT->ProcessLine(".L PRD_DeltaAlphaTInDeltaPTSlices_FSI.cxx++");
	gROOT->ProcessLine("PRD_DeltaAlphaTInDeltaPTSlices_FSI()");

	//gROOT->ProcessLine(".q");		

}
