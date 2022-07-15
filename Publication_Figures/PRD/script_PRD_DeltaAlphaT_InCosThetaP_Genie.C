#include <iostream>
#include <vector>

void script_PRD_DeltaAlphaT_InCosThetaP_Genie() {

	gROOT->ProcessLine(".L PRD_DeltaAlphaTInCosThetaPSlices_Genie.cxx++");
	gROOT->ProcessLine("PRD_DeltaAlphaTInCosThetaPSlices_Genie()");

	gROOT->ProcessLine(".q");		

}
