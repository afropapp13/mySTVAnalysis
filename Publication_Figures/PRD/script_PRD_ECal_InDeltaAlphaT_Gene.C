#include <iostream>
#include <vector>

void script_PRD_ECal_InDeltaAlphaT_Gene() {

	gROOT->ProcessLine(".L PRD_ECalInDeltaAlphaTSlices_Gene.cxx++");
	gROOT->ProcessLine("PRD_ECalInDeltaAlphaTSlices_Gene()");

	gROOT->ProcessLine(".q");		

}
