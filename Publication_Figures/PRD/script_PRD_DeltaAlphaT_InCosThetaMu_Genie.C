#include <iostream>
#include <vector>

void script_PRD_DeltaAlphaT_InCosThetaMu_Genie() {

	gROOT->ProcessLine(".L PRD_DeltaAlphaTInCosThetaMuSlices_Genie.cxx++");
	gROOT->ProcessLine("PRD_DeltaAlphaTInCosThetaMuSlices_Genie()");

	gROOT->ProcessLine(".q");		

}
